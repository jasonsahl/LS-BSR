##
# This is a little framework to make it easier and less error prone to write CGI scripts.
# This defines an object which should be implemented by anyone writing an CGI script.  The object
# is then used to construct the response.
import json
import cgitb
import traceback
from StringIO import StringIO

from twisted.python import reflect

from igs.utils import logging



class CGIPage:
    """
    A class which one implements to create a page that is called via CGI.

    This base class provides some reasonable defaults.  The members are:
    contentType, default = Content-Type: text/html
    headers,     default = {}
    body(self)   - This is a member function which is called and should return a string
    """

    def __init__(self, contentType='Content-Type: text/html', headers=None):
        self.contentType = contentType
        if headers is None:
            self.headers = {}
        else:
            self.headers = headers


    def body(self):
        raise Exception('Please implement me')



def generatePage(cgiPage):
    """
    Takes an instance of CGIPage and generates a page from it,
    sending the proper headers and all that
    """
    cgitb.enable()

    ##
    # A bit evil, I know, but we want all output to go to a logging file
    fout = open('/tmp/webservices.log', 'a')
    logging.OUTSTREAM = fout
    logging.ERRSTREAM = fout
    
    try:
        ##
        # Execute the body first, it may want to add to headers or modify them in soem way as
        # well as contentType
        body = cgiPage.body()
        print cgiPage.contentType
        if cgiPage.headers:
            print '\n'.join([h + ': ' + v for h, v in cgiPage.headers.iteritems()])
        print
        print json.dumps(dict(success=True, data=body))
    except Exception, err:
        print cgiPage.contentType
        print
        stream = StringIO()
        traceback.print_exc(file=stream)
        print json.dumps(dict(success=False, data=dict(stacktrace=stream.getvalue(),
                                                       name=reflect.fullyQualifiedName(reflect.getClass(err)),
                                                       msg=str(err))))
    
