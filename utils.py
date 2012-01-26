"some functions for PyPsy"
import numpy as np
np.set_printoptions(suppress=True)

def dprint(d=None,**kwargs):
    "Pretty-print a dictionary"
    if d is None:
        d = {}
    d.update(kwargs)

    for k,v in d.items():
        if np.squeeze(v).shape == ():
            print k, ' = %.4g'% v
        else:
            print k, ' = \n', v

