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

    
