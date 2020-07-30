"""
Any standard R functions that I am used to using, I could add here.
"""

match = lambda a, b: [b.index(x) if x in b else None for x in a]
"""
returns a vector of the positions of (first) matches of its first argument in its second.
(equivalent of R's match)
"""

def unique(seq):
    """
    uniqify a list without losing order.
    stackoverflow says the hacks in this code make it faster.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
