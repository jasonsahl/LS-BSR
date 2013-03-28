##
# These define useful policies for ensuring a machine is configed properly.  This should probably be off lifted into
# a program/library designed for this, but what I need right now is simple
import os

from igs.utils.logging import errorPrintS
from igs.utils.commands import runSingleProgram, ProgramRunError
from igs.utils.config import configFromMap, configFromStream, configFromEnv, replaceStr


##
# These are default config options, these will be moved to a config file eventually
conf = configFromStream(open('/tmp/machine.conf'), configFromMap({
    'stow': {'package_dir': '/usr/local/stow',
             'base_dir': '/usr/local'},
    'opt': {'package_dir': '/opt/opt-packages',
            'base_dir': '/opt'},
    'config': {'filename': '/tmp/machine.conf'},
    }, configFromEnv()))



##
# Exceptions
class PolicyError(Exception):
    pass


##
# A little helper function
def runSystemEx(cmd):
    """This just ignores all stdout"""
    code = runSingleProgram(cmd, None, errorPrintS)
    if code != 0:
        raise ProgramRunError(cmd, code)



def runInDir(dir, f):
    """Helper function, wraps calling f in a specific directory"""
    curdir = os.getcwd()
    os.chdir(dir)
    f()
    os.chdir(curdir)

def dirExists(dirname):
    """
    Ensure a directory exists, create it if not
    Use fileExists if you want to check for existence but not create
    """
    dirname = replaceStr(dirname, conf)
    if not os.path.exists(dirname):
        try:
            runSystemEx('mkdir -p ' + dirname)
        except:
            raise PolicyError('Could not create directory: ' + dirname)

def fileExists(fname):
    fname = replaceStr(fname, conf)
    if not os.path.exists(fname):
        raise PolicyError('File does not exist: ' + fname)
        
def fileOwner(fname, owner, group=None):
    fname = replaceStr(fname, conf)
    owner = replaceStr(owner, conf)
    who = owner
    if group:
        group = replaceStr(group, conf)
        who += ':' + group

    runSystemEx('chown %s %s' % (who, fname))
        
def dirOwner(dirname, owner, group=None, ignoreError=False):
    """Set owners of a directory, recursively. use fileOwner if you do not want recursive"""
    dirname = replaceStr(dirname, conf)
    owner = replaceStr(owner, conf)
    who = owner
    if group:
        group = replaceStr(group, conf)
        who += ':' + group

    run('chown -R %s %s' % (who, dirname), ignoreError=ignoreError)

def dirPermissions(dirname, perms, ignoreError=False):
    dirname = replaceStr(dirname, conf)
    perms = replaceStr(perms, conf)

    run('chmod -R %s %s' % (perms, dirname), ignoreError=ignoreError)
    
    
def ensurePkg(pkgname):
    """Ensure's a package exists"""
    path = os.path.join(conf('stow.package_dir'), pkgname)
    if not os.path.exists(path):
        raise PolicyError('Package does not exist: ' + path)


def installPkg(pkgname):
    runInDir(conf('stow.package_dir'),
             lambda : runSystemEx('xstow ' + pkgname))

def uninstallPkg(pkgname):
    runInDir(conf('stow.package_dir'),
             lambda : runSystemEx('xstow -D ' + pkgname))    


def pkgFileExists(pkgname, fname):
    """Ensures a file in the package exists"""
    fname = replaceStr(fname, conf)
    fileExists(os.path.join(conf('stow.package_dir'), pkgname, fname))
    

def run(cmd, ignoreError=False):
    """This runs a command, be sure that the command backgrounds, add & if you need to"""
    try:
        runSystemEx(replaceStr(cmd, conf))
    except ProgramRunError:
        if not ignoreError:
            raise
            
    

def executeTemplate(fname):
    """
    Executes a template, at this point this is simple variable substitution in a file
    but it could grow to be more complicated.

    This takes a file that ends in .tmpl and produces a file without the .tmpl with the template
    executed
    """
    fname = replaceStr(fname, conf)
    if not fname.endswith('.tmpl'):
        raise PolicyError('%s does not end in .tmpl' % fname)

    fout = open(fname[:-5], 'w')
    fin = open(fname)

    for line in fin:
        fout.write(replaceStr(line, conf))

    fout.close()
    fin.close()


def executePkgTemplate(pkgname, fname):
    """Just a handy wrapper for executing templates in a specific package"""
    executeTemplate(os.path.join(conf('stow.package_dir'), pkgname, fname))
    

def installOptPkg(pkgname):
    """
    This links an optional package into ${opt.base_dir}/<package-name>

    This cuts off at the last '-' in the pkg name in order to come up with /opt/<package-name>

    If the link already exists, it deletes it and then tries to relink

    For example:
    hadoop-0.20.1 will become ${opt.base_dir}/hadoop
    simple-test-thing-1.2 will become ${opt.base_dir}/simple-test-thing
    """
    outname = '-'.join(pkgname.split('-')[:-1])

    srcName = os.path.join(conf('opt.package_dir'), pkgname)
    dstName = os.path.join(conf('opt.base_dir'), outname)

    if os.path.exists(dstName):
        uninstallOptPkg(pkgname)

    runSystemEx('ln -s %s %s' % (srcName, dstName))
        
def uninstallOptPkg(pkgname):
    outname = '-'.join(pkgname.split('-')[:-1])    
    runSystemEx('rm ' + os.path.join(conf('opt.base_dir'), outname))

def ensureOptPkg(pkgname):
    """Ensure's a package exists"""
    path = os.path.join(conf('opt.package_dir'), pkgname)
    if not os.path.exists(path):
        raise PolicyError('Package does not exist: ' + path)
