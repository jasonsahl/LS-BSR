##
# This is a series of functions for running other programs
# The main piece of code lets you use a generator in order to do sequential thigns
# in a event-loop type function.
# See test and test1 for examples of what this looks like
import sys
import os

import subprocess
from select import select

from igs.utils.logging import logPrint
from igs.utils import functional

##
# How to run a process with subprocess
# pipe = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
#                         shell=True)
# stdin, stdout, stderr = pipe.stdin, pipe.stdout, pipe.stderr


class ProgramRunError(Exception):
    def __init__(self, cmd, code):
        self.cmd = cmd
        self.code = code
        Exception.__init__(self)

    def __str__(self):
        return 'Unable to run program %r with exit code %d' % (self.cmd, self.code)

class ProgramRunner:
    """
    This runs a program.

    The exit status will be in .exitCode
    """

    def __init__(self, cmd, stdoutf, stderrf, addEnv=None, env=None, log=False):
        """
        addEnv takes the contents of addEnv and adds them to the current environment.  env only passes what is specified
        as the environment
        """
        self.cmd = cmd
        self.stdoutf = stdoutf
        self.stderrf = stderrf
        self.addEnv = addEnv
        self.env = env
        self.log = log
        self.exitCode = None

    def __call__(self):
        """
        This returns:
        (onComplete, [(stream1, func1), .. (streamn, funcn)])

        Where onComplete is any cleanup that needs to happen once all the streams are consumed
        and stream1 is a stream and func1 is the function to call upon data coming in for that stream
        """
        if self.log:
            logPrint(self.cmd)


        env = self.env
        if self.addEnv and not self.env:
            # Copy the current environment because we'll be modifying it
            env = functional.updateDict(dict(os.environ), self.addEnv)
        elif self.addEnv:
            env = functional.updateDict(dict(self.env), self.addEnv)
            
        pipe = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                shell=True, env=env)
        self.pipe = pipe
                                
        return (self.onComplete, [(pipe.stdout, self.stdoutf), (pipe.stderr, self.stderrf)])

    def onComplete(self):
        self.exitCode = self.pipe.wait()
        
        self.pipe.stdout.close()
        self.pipe.stderr.close()
        #self.pipe.stdin.close()



def runProgramRunner(pr):
    """Runs a ProgramRunner, blocking until it is finished returning the exit code"""
    def _():
        yield pr

    runCommandGens([_()])

    return pr.exitCode

def runProgramRunnerEx(pr):
    code = runProgramRunner(pr)

    if code != 0:
        raise ProgramRunError(pr.cmd, code)

def runSingleProgram(cmd, stdoutf, stderrf, addEnv=None, env=None, log=False):
    """Gives you control over where the stream data goes"""
    pr = ProgramRunner(cmd, stdoutf, stderrf, addEnv=addEnv, env=env, log=log)

    def _():
        yield pr

    runCommandGens([_()])

    return pr.exitCode

def runSingleProgramEx(cmd, stdoutf, stderrf, addEnv=None, env=None, log=False):
    """Gives you control over where the stream data goes"""
    pr = ProgramRunner(cmd, stdoutf, stderrf, addEnv=addEnv, env=env, log=log)

    def _():
        yield pr

    runCommandGens([_()])

    if pr.exitCode != 0:
        raise ProgramRunError(pr.cmd, pr.exitCode)
    
    return pr.exitCode


def runSystem(cmd, addEnv=None, env=None, log=False):
    """Simpler wrapper, looks more like os.system"""
    return runSingleProgram(cmd, sys.stdout.write, sys.stderr.write, addEnv=addEnv, env=env, log=log)

def runSystemEx(cmd, addEnv=None, env=None, log=False):
    """Wrapper around runSystem, throws an exception with error code on non zero return"""
    code = runSystem(cmd, addEnv=addEnv, env=env, log=log)
    if code != 0:
        raise ProgramRunError(cmd, code)
        
        
def getStreams(state):
    (_gen, (_oncomplete, streams)) = state

    return streams

def getOnComplete(state):
    (_gen, (onComplete, _streams)) = state

    return onComplete

def ctorGenerators(gens):
    """Construct the generators"""

    return [iter(g) for g in gens]


def activeOutputStreams(states):
    """Returns a list of all the active output streams and their
    index in 'states'"""
    
    res = []
    
    for idx, s in enumerate(states):
        if s:
            _, v = s
            if v:
                _, streams = v
                res.extend([(x, idx) for x in streams.keys()])

    return res

def nextIteration(states):
    """
    This actually modifies states, it just returns
    the modified states to be more functional.  THIS MAY CHANGE TO NOT MODIFY STATES
    BUT PRODUCE A NEW VALUE
    """
    
    for idx, s in enumerate(states):
        if s:
            g, v = s
            if not v:
                try:
                    runner = g.next()
                    onComplete, streams = runner()
                    states[idx] = (g, (onComplete, dict(streams)))
                except StopIteration:
                    # Turn this guy off if we are done
                    states[idx] = None

    return states



def callStreamF(stream, data, state):
    streams = getStreams(state)

    if streams[stream]:
        streams[stream](data)

    
def removeStream(stream, state):
    streams = getStreams(state)
    del streams[stream]

    ##
    # Really returning this beacuse callers want to know if there
    # are any more streams in there. might as well giev them
    # all the streams incase they want to do anything with it
    return streams
    
        
def runCommandGens(generators):
    """
    This takes a list of generators and runs through each of them
    in parallel running the commands.

    Each iteration of a generator must return something that is callable and returns a tuple that looks like:
    (oncomplete, [(stream1, function), (stream2, function), ... (streamn, function)])

    Where 'oncomplete' is a function that gets called when all streams have been exhausted
    and each stream has an associated function that gets called with the output
    """

    ##
    # contain objects being worked on from generator
    states = [(g, None) for g in ctorGenerators(generators)]

    ##
    # start initial commends:
    states = nextIteration(states)

    outputStreams = dict(activeOutputStreams(states))

    while outputStreams:
        input, _output, _error = select(outputStreams.keys(), [], [])

        iterateAndBuild = False
        for s in input:
            line = s.readline()

            if line:
                callStreamF(s, line, states[outputStreams[s]])
            else:
                iterateAndBuild = True
                            
                ##
                # removeStream returns True if this was the last stream to be removed
                # which means we need to try to start up the next iteration
                # and recreate outputStreams
                if not removeStream(s, states[outputStreams[s]]):
                    state = states[outputStreams[s]]
                    f = getOnComplete(state)
                    ##
                    # Finished with this one, call onComplete
                    f()
                    # Set the second portion to None, nextIteration uses that to know
                    # when to get the next iteration
                    states[outputStreams[s]] = (state[0], None)

        # If any of our streams are completely done, lets get new ones and rebuild
        if iterateAndBuild:
            states = nextIteration(states)
            outputStreams = dict(activeOutputStreams(states))
                


def test():
    import random
    
    def _(id):
        print id
        yield ProgramRunner('date', sys.stdout.write, sys.stderr.write, log=True)
        print id
        yield ProgramRunner('sleep %d' % random.randrange(0, 10), sys.stdout.write, sys.stderr.write, log=True)
        print id
        yield ProgramRunner('date', sys.stdout.write, sys.stderr.write, log=True)
        print id


    runCommandGens([_(v) for v in range(100)])
    

def test1():

    def f1():
        yield ProgramRunner('date', sys.stdout.write, sys.stderr.write, log=True)
        yield ProgramRunner('echo what is up', sys.stdout.write, sys.stderr.write, log=True)
        yield ProgramRunner('sleep 3', sys.stdout.write, sys.stderr.write, log=True)        
        yield ProgramRunner('echo that is cool', sys.stdout.write, sys.stderr.write, log=True)

    def f2():
        yield ProgramRunner('sleep 2', sys.stdout.write, sys.stderr.write, log=True)
        yield ProgramRunner('echo this is another thing running concurrently', sys.stdout.write, sys.stderr.write, log=True)


    runCommandGens([f1(), f2()])
    
