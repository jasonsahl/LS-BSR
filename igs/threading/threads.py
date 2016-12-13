##
# Some useful functions for threads
import threading

from igs.utils import functional as f

from igs.threading.channels import Channel


def runThread(func, *args, **kwargs):
    th = threading.Thread(target=func, args=args, kwargs=kwargs)
    th.start()
    return th

def runThreadWithChannel(func):
    """
    This creates a thread running the passed function.  The only argument to the function
    is the channel that will be used to communicate with the thread.
    """
    ch = Channel()
    th = runThread(func, ch)
    return f.Record(thread=th, channel=ch)

def mp_shell(func, params, numProc):
    from multiprocessing import Pool
    from functools import partial
    pool = Pool(processes=int(numProc))
    new_func = partial(func,params)
    pool.map(new_func,params)
    pool.terminate()
