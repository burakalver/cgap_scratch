'''
Given a search result from snovault, I would like flatten the results out.
e.g. take a list with 2 elements
[
  {
    'a':1,
    'b1':{'c':1},
  },
  {'b1':{'c':2},
    'b2':[{'d1':3,'d2':4},{'d2':5}]
  }
]

let's study the keys and values for this
{
   a:1,
   b1.c:[1,2]
   b2.d1:3
   b2.d2:[4,5]
}

'''


def nested_keys(parent0, depth_max=-1):
    '''
    parent0 is a list or dict with possible deep embedddings of lists and dicts.
    Returns list of keys, in the flattened format. e.g. [a, b1.c, b2.d1, b2.d2]
    '''
    allkeys = set()
    parentkey = []
    depth = 0
    nested_keys_routine(parent0, allkeys, parentkey, depth, depth_max)
    return allkeys


def nested_keys_routine(parent, allkeys, parentkey, depth, depth_max):
    ''''function to be called recursively for nested_keys '''
    if depth <= depth_max or depth_max == -1:
        if isinstance(parent, list):
            for val in parent:
                nested_keys_routine(val, allkeys, parentkey, depth+1, depth_max)
        if isinstance(parent, dict):
            for key, val in sorted(parent.items()):
                currentkey = parentkey.copy()
                currentkey.append(key)
                if not isinstance(val, dict) and not isinstance(val, list):
                    allkeys.add(".".join(currentkey))
                nested_keys_routine(val, allkeys, currentkey, depth+1, depth_max)


def get_field_by_nested_key(parent0, keystr):
    '''
    parent is a list or dict with possible deep embedddings of lists and dicts.
    the keystr is "."-separated list of keys down the embedding.
    if parent=[{'a':1},{'b1':{'c':2},'b2':[{'d1':3,'d2':4},{'d2':5}]}]
    and keystr = b2.d2
    returns [4,5]
    '''
    keys = keystr.split(".")
    values = []
    get_field_by_nested_key_routine(parent0, keys, values)
    return values


def get_field_by_nested_key_routine(parent, keys, values):
    '''function to be called recursively for get_field_by_nested_key '''
    
    if len(keys) == 0:
        return
    if isinstance(parent, list):
        for val in parent:
            get_field_by_nested_key_routine(val, keys, values)
    else:
        remainingkeys = keys.copy()
        curkey = remainingkeys.pop(0)
        val = parent.get(curkey, None)
        if isinstance(val, dict):
            get_field_by_nested_key_routine(val, remainingkeys, values)
        elif isinstance(val, list):
            if (not isinstance(val[0], (list, dict))) and len(remainingkeys)==0:
                values.extend(val)
            else:
                get_field_by_nested_key_routine(val, remainingkeys, values)
        else:
            if val is not None:
                values.append(val)
