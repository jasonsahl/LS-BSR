##
# Channels provide a means of communicatin between threads
from Queue import Queue

class Channel:
    """
    A channel is a unidirectional form of communication between threads.  
    
    A channel allows for the following actions:

    send - Send an object over the channel
    sendError - Sends an error, the object should be an exception.  It will get raised on the otherside when they do a receive
    receive - Receive an object from a channel, this blocks unless a timeout is set
    NOTE - One thread should only be send'ing and the other thread only receive'ing.  Channels unidirectional

    sendWithChannel - Send an object and create a receive channel, and return it.  The object will be sent as a tuple (object, Channel).
                      This is useful if you want to send a task to perform and get a result back, the channel is almost like a 'future'
                      but not quite. There is no sendErrorWithChannel.


    """
    
    def __init__(self):
        self.queue = Queue()


    def send(self, obj):
        """Send 'obj' through the channel"""

        self.queue.put_nowait((True, obj))


    def sendError(self, err):
        """Send an error"""
        self.queue.put_nowait((False, err))
        
    def sendWithChannel(self, obj):
        """Send 'obj' as well as a new channel through this channel and return the new channel"""
        ch = Channel()
        self.send((obj, ch))
        return ch

        
    def receive(self, timeout=None):
        """
        Receive an object from the channel.  Blocks unless a timeout is specified.  Give a timeout of 0 to poll
        """

        block = True
        
        ok, item = self.queue.get(block, timeout)
        self.queue.task_done()
        if not ok:
            raise item
        return item
