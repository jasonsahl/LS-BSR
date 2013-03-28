##
# Little tools to make working with the CLI easier
import optparse

from igs.utils.config import configFromStream, configFromMap, configFromEnv, replaceStr

from igs.utils.functional import applyIfCallable

##
# These are the different control types we can take as input.  For instance making a list, counting, string
# string is the default.
# BINARY is True for backwards compatibility
BINARY = True
STRING = 'string'
LIST = 'list'
COUNT = 'count'


class DeprecatedOptionError(Exception):
    pass

class MissingOptionError(Exception):
    pass

class InvalidOptionError(Exception):
    pass

class CLIError(Exception):
    def __init__(self, option, original):
        self.msg = str(original)
        self.option = option

    def __str__(self):
        return 'Error handling option: %s, failed with message: %s' % (self.option, self.msg)



def applyOption(val, option, conf):
    try:
        #
        # We want to apply any replacements on the options
        # The question is if baseConf is really the config file
        # we should be applying these from...
        v = option[4](val)
    
        #
        # Now v might still be a function at this point because we want to give
        # them the ability to access baseConf is they need it to do more replacements
        v = applyIfCallable(v, conf)
        try:
            return replaceStr(v, conf)
        except TypeError:
            return v
    except MissingOptionError, err:
        raise CLIError(option[0], err)
    
    
def buildConfigN(options, args=None, usage=None, baseConf=None, putInGeneral=True):
    """
    This builds a config from options.  Options is a list of tuples that looks like:

    (name, short, long, help, func, [bool])
    Where
    name - Name of the option, this is what it will become in the config file
    short - Short option - needs to start with -
    long - Long option - needs to start with --
    help - Help to be given to a user in --help output
    func - Function to be applied to the value
    bool - This is not required, set to True if the option is simply a boolean, all other datatypes can be verified via 'func'

    This will implicitly check if a 'conf' option exists, and if so load te conf file as a base for these config options.

    All options are put into the 'general' section.

    This returns a tuple
    (conf, args)
    where args is whatever is left over from parsing

    This also implicitly loads the current environment into the env section

    If, when evaluated, 'func' returns a function, it is called with the baseConf.  This is to allow more complex replacements to
    happen.
    """
    def _iterBool(v):
        """
        Adds the non erquired bool field with a default of STRING if
        it is not present
        """
        for l in v:
            if len(l) == 6:
                yield l
            else:
                yield tuple(list(l) + [STRING])
                
    
    parser = optparse.OptionParser(usage=usage)

    ##
    # keep track of the function to apply to conf
    confFunc = None

    ##
    # The order of the options below
    # var name, short option, long option, help message, function to apply, binary option    
    for n, s, l, h, f, b in _iterBool(options):
        ##
        # We could have a function we want to apply to the conf variable.  We want to store it
        # so when we use it in the next block we don't have to loop over options looking for it again
        # This is a minor optimization and probably not even necessary...
        if n == 'conf':
            confFunc = f
            
        if b == BINARY:
            parser.add_option(s, l, dest=n, help=h, action='store_true')
        elif b == LIST:
            parser.add_option(s, l, dest=n, help=h, action='append')
        elif b == COUNT:
            parser.add_option(s, l, dest=n, help=h, action='count')
        elif b == STRING:
            parser.add_option(s, l, dest=n, help=h)
        else:
            raise Exception('Unknown option type: ' + repr(b))

    ops, args = parser.parse_args(args=args)

    if baseConf is None:
        baseConf = configFromEnv()
    if hasattr(ops, 'conf'):
        baseConf = configFromStream(open(replaceStr(confFunc(ops.conf), baseConf)), baseConf)

    vals = {}

    ##
    # The order of the options below
    # var name, short option, long option, help message, function to apply, binary option
    for o in _iterBool(options):
        n, _s, l, _h, f, _b = o
        try:
            vals[n] = applyOption(getattr(ops, n), o, baseConf)
        except Exception, err:
            raise CLIError(l, err)
            

    if putInGeneral:
        vals = {'general': vals}
    return (configFromMap(vals, baseConf), args)



##
# These are various functions to make building and verifying data easier
def notNone(v):
    """
    Throws MissingOptionError if v is None, otherwise returns v
    """
    if v is None:
        raise MissingOptionError('Must provide a value for option')

    return v


def defaultIfNone(d):
    """
    Returns a function that returns the value 'd' if the passed
    value to the function is None
    """
    def _(v):
        if v is None:
            return d
        else:
            return v

    return _

def restrictValues(values):
    def _(v):
        if v not in values:
            raise InvalidOptionError('Value must be one of: %s' % ', '.join([str(x) for x in values]))
        return v

    return _

def notFalse(v):
    """
    Throws an exception if option is blank
    """
    if not v:
        raise MissingOptionError('Must provide a value for option')

    return v

# This is just a more descriptive name for strings than notFalse
notBlank = notFalse

def deprecated(msg):
    def _(x):
        if x is not None:
            raise DeprecatedOptionError(msg)
        else:
            return x

    return _

def composeCLI(*funcs):
    """
    This function is like compose except inbetween each function
    it does a replaceStr from a config on the intermediate values
    if it is a string.  Usage:
    composeCLI(f, g)(x)(conf)
    """
    funcs = list(funcs)
    funcs.reverse()
    def v(x):
        def c(conf):
            val = x
            for f in funcs:
                val = f(val)
                try:
                    val = replaceStr(val, conf)
                except TypeError:
                    pass

            val = applyIfCallable(val, conf)
            try:
                return replaceStr(val, conf)
            except TypeError:
                return val
        return c
    return v
