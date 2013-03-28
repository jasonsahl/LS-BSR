##
# This is for working with ini style configs.  This offers a few tricks though:
#
# 1) You can reference variables in this by ${vname}
#    For example you can have:
#    foo = bar
#    bar = baz
#    foobar = ${foo}${bar}
#    And foobar will be barbaz
# 2) The name of a section is really just a short hand for writing fully dot separated names out. For example:
#    [general]
#    bar = baz
#    foo = zoom
#
#    This is really equivalent to:
#    []
#    general.bar = baz
#    general.foo = zoom
#
#    Note the [] at the top says we went the root section.  This means you can access variables as: ${general.bar} if you want to access it in another section.
#    Names in a particular section are searched in that section and then search the root name.  For example:
#    []
#    foo = bar
#    [general]
#    baz = zoom
#    t1 = ${baz}
#    t2 = ${foo}
#
#    t1 is zoom, t2 is bar
# 3) You can also add to sections defined elsewhere in the file (although not suggested).  For example:
#    [general]
#    foo = bar
#    [other]
#    test = testing
#    [general]
#    bar = baz
#
# 4) An include facility is provided.  This is through 2 directives: -include and -include_in_section.  The semantics are as follows:
#    -include means to read in the provided file and to include its definitions irrespective of any context in the current file.
#    For example, if you have f1.conf defined as:
#    [foo]
#    bar=baz
#
#    Then f2.conf defined as:
#    [boom]
#    -include f1.conf
#    zoom=zoop
#
#    This would be equivalent to the following file
#    [boom]
#    [foo]
#    bar=baz
#    [boom]
#    zoom=zoop
#
#    -include_in_section adds the contectual information of the file-including the other.  If f1.conf is defined the same as
#    above, but f2.conf is defined as:
#    [boom]
#    -include_in_section f1.conf
#    zoom=zoop
#
#    That is equivalent to the following file:
#    [boom]
#    foo.bar=baz
#    zoom=zoop
#
#    NOTE: -include and -include_in_section are *NOT* like #include in C/C++.  This is not a basic copy and paste mechanism.  Instead they
#    have contextual knowledge of where the processing currently is in the config file and apply that knowledge during the include.
import os
import re

##
# Regex to find variables in a string so we can replace on them
VARIABLE_RE = re.compile(r'\${([A-Za-z0-9_\.]+)}')

class ConfigParseError(Exception):
    pass

class NoKeyFoundError(Exception):
    def __init__(self, k):
        self.key = k
        Exception.__init__(self)

    def __str__(self):
        return self.key

###
# New config style, more flexible
class Config:
    """
    This represents a config, it's basically a linked list, each
    config pointing to its parent and parent being None if there is none

    You'll notice there are no methods to modify this, Configs are meant to be
    immutable
    """

    def __init__(self, m, base, lazy=False, depth=None):
        self.conf = flattenMap(m, depth=depth)
        self.base = base
        self.lazy = lazy


        ##
        # If lazy is set, don't validate
        if not lazy:
            for k in self.conf.keys():
                v = self.conf[k]
                newV = replaceVariables(k, v, self)
                #if newV != v:
                #    self.conf[k] = newV


    def __call__(self, key, **kwargs):
        """
        We want this to match the current API which looks like function calls
        """
        return self.get(key, **kwargs)

    def __contains__(self, key):
        return key in self.keys()

    def get(self, key, **kwargs):
        return replaceVariables(key,
                                self.get_raw(key, **kwargs),
                                self.get_raw)
    
    def get_raw(self, key, **kwargs):
        """
        If the key does not exist in this config and there is no base then raise
        NoKeyFoundError

        If the key does not exist and 'default' is in kwargs and there is no base
        then return 'default'

        If the key does not exist and there is a base then ask the base

        Otherwise, return the key
        """
        if key not in self.conf and not self.base and 'default' not in kwargs:
            raise NoKeyFoundError(key)
        elif key not in self.conf and not self.base and 'default' in kwargs:
            return kwargs['default']
        elif key not in self.conf:
            return self.base.get_raw(key, **kwargs)
        else:
            return self.conf[key]


        
    def keys(self):
        """
        Return a set of all the keys
        """
        s = set(self.conf.keys())
        if self.base:
            s = s | self.base.keys()

        return s
    

def configListFromStream(stream):
    """
    Constructs a list of [(key, value)] from an ini-like file
    """
    ##
    # We are just going to read the config file into a map and use configFromMap
    cfg = []
    ##
    # Default section
    section = ''
    for line in stream:
        line = line.lstrip()
        if not line or line[0] == '#':
            # A comment or empty line
            continue
        elif line.startswith('-include '):
            _, path = line.split(' ', 1)
            path = path.strip()
            cfg.extend(configListFromStream(open(path)))
            continue
        elif line.startswith('-include_in_section '):
            _, path = line.split(' ', 1)
            path = path.strip()
            if section:
                cfg.extend([(section + '.' + k, v) for k, v in configListFromStream(open(path))])
            continue
        
        if line[0] == '[':
            line = line.rstrip()
            if line[-1] == ']':
                section = line[1:-1]
                continue
            else:
                raise ConfigParseError("""Error parsing line: %r""" % line)
        else:
            key, value = line.split('=', 1)
            ##
            # Remove the trailing '\n'
            if section:
                key = section + '.' + key
            cfg.append((key, value.rstrip('\n')))

    return cfg
    
    
def configFromStream(stream, base=None, lazy=False, depth=None):
    """
    Constructs a config function from a stream.

    base is used if you want to make this on top of base (for example this references variables in base.

    The returned value will be a composite of base with stream."""

    return configFromMap(dict(configListFromStream(stream)), base, lazy, depth)


def flattenMap(map, depth):
    res = {}
    for k, v in map.iteritems():
        if k:
            n = k + '.'
        else:
            n = ''
        #
        # Cheap way to see if v is a map
        # We only want to do this if depth is true (a positive int)
        # or depth is None (infinit depth)
        if hasattr(v, 'keys') and (depth > 0 or depth is None):
            if depth: depth = depth - 1
            for kp, vp in flattenMap(v, depth).iteritems():
                res[n + kp] = vp
        else:
            res[k] = v

    return res

def replaceVariables(k, value, lookup):
    """
    Takes value and determines if any variables need replacing in it and replaces
    them with whatever is in 'lookup'.  This will blow out the stack if
    you have a recursive variables.

    If a variable cannot be found a NoKeyFoundError is thrown.
    """
    ##
    # We want the section for this variable which means if the key is
    # a.b.c.d the section is a.b.c
    # If the key is empty, the section is empty
    if '.' in k:
        ##
        # Add trailing . so we don't have to later
        section = '.'.join(k.split('.')[:-1]) + '.'
    else:
        section = ''

    ##
    # Check to make sure this is a string:
    if hasattr(value, 'replace'):
        vars = VARIABLE_RE.findall(value)
        if vars:
            for var in vars:
                if '.' not in var:
                    secVar = section + var
                else:
                    secVar = var

                try:
                    value = value.replace('${%s}' % var, lookup(secVar))
                except NoKeyFoundError:
                    value = value.replace('${%s}' % var, lookup(var))

            return replaceVariables(secVar, value, lookup)
        else:
            return value
    else:
        return value


def replaceStr(str, lookup):
    """Takes a string and replaces all variables in it with values from config"""
    vars = VARIABLE_RE.findall(str)
    if vars:
        for var in vars:
            str = str.replace('${%s}' % var, lookup(var))

    return str

def configFromMap(map, base=None, lazy=False, depth=None):
    return Config(map, base, lazy, depth)


def configFromDict(d, base=None, lazy=False, depth=None):
    return configFromMap(d, base, lazy, depth)


def configFromConfig(c, base=None, lazy=False, depth=None):
    return configFromMap(configToDict(c), base, lazy, depth)

def configFromEnv(base=None, lazy=False, depth=None):
    """
    This constructs a config from the environment variables and puts them in the [env] section.
    Case is preserved, so it would be ${env.HOME}.  Variables are also replaced.
    """
    return configFromMap({'env': os.environ}, base, lazy, depth)


def configToDict(config):
    if config.lazy:
        return dict([(k, config.get_raw(k)) for k in config.keys()])
    else:
        return dict([(k, config.get(k)) for k in config.keys()])

def test():
    c = configFromMap({
        'general': {'foo': '${bar}',
                    'bar': 'baz'},
        'test': {'foo': '${general.foo}',
                 'bar': '${zoom}',
                 'zoom': 'bye'},
        })


    assert 'baz' == c('general.foo')
    assert 'baz' == c('general.bar')
    assert 'baz' == c('test.foo')
    assert 'bye' == c('test.bar')
    assert 'bye' == c('test.zoom')


    c1 = configFromMap({
        'general': {'zoom': '${bar}'},
        'testing': {'zoom': '${test.zoom}'}
        }, c)

    assert 'baz' == c1('general.zoom')
    assert 'bye' == c1('testing.zoom')
    
    
