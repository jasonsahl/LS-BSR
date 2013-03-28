##
# Functions for creating and reading a request
import cgi
import json
import urllib
import httplib
import socket

from igs.utils.errors import TryError

def performQueryNoParse(host, url, var, timeout=30, debug=False):
    def _performQuery():
        params = urllib.urlencode({'request': json.dumps(var)})
        conn = httplib.HTTPConnection(host, timeout=timeout)
        conn.connect()
        if debug:
            conn.set_debuglevel(3)
        ##
        # Incredibly cheap hack
        conn.sock.settimeout(60)
        conn.request('POST', url, params, headers={'Content-Type': 'application/x-www-form-urlencoded'})
        return conn.getresponse().read()

    count = 4
    while True:
        try:
            return _performQuery()
        except socket.timeout:
            count -= 1
            if count <= 0:
                raise
    

def performQuery(host, url, var, timeout=30, debug=False):
    """
    params is a dict on of values to pass to server
    """
    rawData = performQueryNoParse(host, url, var, timeout=timeout, debug=debug)
    try:
        result = json.loads(rawData)
        data = result['data']
        if not result['success']:
            raise TryError('Query failed: ' + data['msg'], result)
        return data
    except TryError:
        raise
    except Exception:
        raise ValueError('Unknown data: ' + str(rawData))

def readQuery():
    form = cgi.FieldStorage()
    return json.loads(form['request'].value)
    
