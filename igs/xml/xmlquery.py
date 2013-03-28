
def prop(name, value):
    """Simple implementation of checking a property, will need to add
    better pattern checking
    """
    def _(elem):
        if hasattr(elem, 'hasAttribute'):
            return elem.hasAttribute(name) and elem.getAttribute(name) == value
        else:
            return False
    return _

def name(name):
    """Checks if the element has a particular name"""
    def _(elem):
        if hasattr(elem, 'tagName'):
            return elem.tagName == name
        else:
            return False
    return _

def execQuery(query, doc):
    """Executs a query on a document"""
    ##
    # Determine if dealing with a Document or Element
    if hasattr(doc, 'doctype'):
        return _execQuery(query, doc.childNodes)
    else:
        return _execQuery(query, [doc])

def _execQuery(query, elems):
    nodes = list(elems)
    noncall = []
    for q in query:
        try:
            nodes = filter(q, nodes)
        except TypeError:
            noncall.extend(q)

    res = []
    for n in nodes:
        res.extend(n.childNodes)
    if noncall:
        return _execQuery(noncall, res)
    else:
        return res


def test():
    from xml.dom import minidom

    doc = minidom.parseString("""
    <slideshow>
    <title>Demo slideshow</title>
    <slide><title>Slide title</title>
    <point>This is a demo</point>
    <point>Of a program for processing slides</point>
    </slide>
    
    <slide><title>Another demo slide</title>
    <point>It is important</point>
    <point>To have more than</point>
    <point>one slide</point>
    </slide>
    </slideshow>
    """)

    query = [name('slideshow'),
             [name('slide')]
             ]

    return execQuery(query, doc)
    
