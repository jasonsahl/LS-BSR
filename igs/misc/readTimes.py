## A little utility to parse run times out of pipeline XML files

import sys
import time
import math

from xml.dom import minidom

from igs.xml.xmlquery import execQuery, name


def chunk(i, ch):
    chunk = []
    for v in i:
        chunk.append(v)
        if len(chunk) == ch:
            yield chunk
            chunk = []

    if chunk:
        yield chunk

def parseTime(t):
    """
    Want to parse something that looks like:
    2010-05-25T08:14:00.357Z
    """
    return time.mktime(time.strptime(t.split('.')[0], '%Y-%m-%dT%H:%M:%S'))


def parseFile(f):
    doc = minidom.parse(f)
    query = [name('commandSetRoot'),
             [name('commandSet'),
              [name('command')]]]

    res = execQuery(query, doc)
    
    times = [r.childNodes[0].data for r in res if r.localName in ['startTime', 'endTime']]

    return [parseTime(end) - parseTime(start) for start, end in chunk(times, 2)]
    

def getMax(td):
    every = []
    for f, v in td.iteritems():
        every.extend([(v1, f) for v1 in v])

    (v, f) = max(every)

    return (f, v)

def getMin(td):
    every = []
    for f, v in td.iteritems():
        every.extend([(v1, f) for v1 in v])

    (v, f) = min(every)

    return (f, v)

def getAvg(td):
    every = []
    for _, v in td.iteritems():
        every.extend(v)

    return sum(every)/len(every)

def getStdDev(td):
    every = []
    for _, v in td.iteritems():
        every.extend(v)
        
    avg = getAvg(td)

    return math.sqrt(sum([(e - avg)**2 for e in every])/len(every))

def main():
    timingData = {}
    for f in sys.stdin:
        f = f.strip()
        ##
        # Only take stuff greater than no time
        timingData[f] = [t for t in parseFile(f) if t > 3.0]


    maxFile = getMax(timingData)
    minFile = getMin(timingData)
    avgTime = getAvg(timingData)
    stdDev = getStdDev(timingData)

    print 'maxFile: %s %f' % maxFile
    print 'minFile: %s %f' % minFile
    print 'avg: %f' % avgTime
    print 'stdDev: %f' % stdDev
        

if __name__ == '__main__':
    main()
