##
# Useful functions for doing things with ssh
from igs.utils.logging import errorPrint, errorPrintS
from igs.utils.commands import ProgramRunError, ProgramRunner, runProgramRunner
from igs.utils import commands

def quoteStr(s):
    """
    Simple implementation right now
    """
    return "'" + s + "'"

def runSystemSSHA(host, cmd, stdoutf, stderrf, user=None, options=None, log=False):
    """
    Asynchornous function for running a command through ssh
    """
    command = ['ssh']
    if user:
        host = user + '@' + host

    if options:
        command.append(options)

    command.append(host)

    command.append(quoteStr(cmd))

    return ProgramRunner(' '.join(command), stdoutf, stderrf, log=log)

def runSystemSSH(host, cmd, stdoutf, stderrf, user=None, options=None, log=False):
    """
    Blocking version of runSystemSSHA

    Returns the exit code
    """
    return runProgramRunner(runSystemSSHA(host, cmd, stdoutf, stderrf, user, options, log=log))


def runSystemSSHEx(host, cmd, stdoutf, stderrf, user=None, options=None, log=False):
    """
    Blocking version of runSystemSSHA, throws an exception on failure though
    """
    pr = runSystemSSHA(host, cmd, stdoutf, stderrf, user, options, log=log)
    exitCode = runProgramRunner(pr)
    if exitCode != 0:
        raise ProgramRunError(pr.cmd, exitCode)


def scpFromA(host, src, dst, user=None, options=None, log=False):
    """
    Asynchronous function to copy scp a file from a host
    """
    command = ['scp']
    if user:
        host = user + '@' + host

    if options:
        command.append(options)

    command.append(host + ':' + src)
    command.append(dst)

    return ProgramRunner(' '.join(command), None, errorPrint, log=log)

def scpToA(host, src, dst, user=None, options=None, log=False):
    """
    Asynchronous function to copy scp a file to a host
    """
    command = ['scp']
    if user:
        host = user + '@' + host

    if options:
        command.append(options)

    command.append(src)
    command.append(host + ':' + dst)

    return ProgramRunner(' '.join(command), None, errorPrint, log=log)

def scpFrom(host, src, dst, user=None, options=None, log=False):
    """
    Blocking function to copy scp a file from a host

    Returns the exit code
    """
    return runProgramRunner(scpFromA(host, src, dst, user, options, log=log))

def scpFromEx(host, src, dst, user=None, options=None, log=False):
    """
    Blocking function to copy scp a file from a host

    Throws an exception if it fails
    """
    pr = scpFromA(host, src, dst, user, options, log=log)
    exitCode = runProgramRunner(pr)
    if exitCode != 0:
        raise ProgramRunError(pr.cmd, exitCode)

def scpTo(host, src, dst, user=None, options=None, log=False):
    """
    Blocking function to copy scp a file to a host

    Returns the exit code
    """
    return runProgramRunner(scpToA(host, src, dst, user, options, log=log))

def scpToEx(host, src, dst, user=None, options=None, log=False):
    """
    Blocking function to copy scp a file to a host

    Throws an exception if it fails
    """
    pr = scpToA(host, src, dst, user, options, log=log)
    exitCode = runProgramRunner(pr)
    if exitCode != 0:
        raise ProgramRunError(pr.cmd, exitCode)    

def rsyncFromEx(host, src, dst, rsyncOptions=None, user=None, log=False):
    cmd = ['rsync']
    if rsyncOptions:
        cmd.append(rsyncOptions)

    if user:
        host = user + '@' + host

    cmd.extend([host + ':' + src, dst])

    return commands.runSingleProgramEx(' '.join(cmd), stdoutf=None, stderrf=errorPrintS, log=log)

def rsyncFrom(host, src, dst, rsyncOptions=None, user=None, log=False):
    try:
        rsyncFromEx(host, src, dst, rsyncOptions, user, log)
        return 0
    except ProgarmRunError, err:
        return err.exitCode

           
           
