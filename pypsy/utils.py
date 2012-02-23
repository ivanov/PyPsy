"some functions for PyPsy"
import numpy as np
import scipy.special as special
np.set_printoptions(suppress=True)

def dprint(d=None,**kwargs):
    """Pretty-print a dictionary.

    You can also pass additional keyword arguments which will get placed in
    the dictionary that's printed.

    Examples
    --------
    >>> mat2d = np.atleast_2d(np.arange(0,1,.2)).T
    >>> array1d = np.arange(0,1,.3333) - .00000001
    >>> d=dict(a_really_long_name=mat2d, a=1.3, b=.33, x=.2, bar=array1d)
    >>> print d # random: dictionaries have no guaranteed order
    {'a': 1.3, 'x': 0.20000000000000001, 'a_really_long_name': array([[ 0. ],
           [ 0.2],
           [ 0.4],
           [ 0.6],
           [ 0.8]]), 'b': 0.33000000000000002, 'bar': array([-0.00000001,  0.33329999,  0.66659999,  0.99989999])}
    >>> dprint(d)
    a = 1.3
    a_really_long_name = [[ 0. ]
                          [ 0.2]
                          [ 0.4]
                          [ 0.6]
                          [ 0.8]]
    b = 0.33
    bar = [-0.00000001  0.33329999  0.66659999  0.99989999]
    x = 0.2

    >>> dprint(d, another_var='hello')
    a = 1.3
    a_really_long_name = [[ 0. ]
                          [ 0.2]
                          [ 0.4]
                          [ 0.6]
                          [ 0.8]]
    another_var = "hello"
    b = 0.33
    bar = [-0.00000001  0.33329999  0.66659999  0.99989999]
    x = 0.2
    """
    if d is None:
        d = {}
    d.update(kwargs)

    for k,v in sorted(d.items()):
        if isinstance(v, basestring):
            print k, '= "%s"'% v
        elif np.squeeze(v).shape == ():
            print k, '= %.4g'% v
        else:
            line = k + ' ='
            vrep = str(v)
            # pad every line with space to match the name
            vlines = ("\n" + " " *len(line) + " ").join(vrep.splitlines())
            print line, vlines

def ztop(z):
    p = 0.5 + 0.5*special.erf(z/np.sqrt(2))
    return p

def ptoz(p):
    p = np.asanyarray(p)
    z = np.sqrt(2)*special.erfinv(p*2-1);
    return z

def afc2_p_to_dprime(ph, pf=None):
    if pf is None:
        pf = ph[:,1]
        ph = ph[:,0]

    zH=ptoz( ph )
    zF=ptoz( pf )

    dP=(zH-zF)/np.sqrt(2)
    C=-(zH+zF)/np.sqrt(2)
    lnB=(zF**2-zH**2)/2
    pC=(ph+(1.0-pf))/2

    return dP, C, lnB, pC

def afc2_dprime_to_p(dp, c):
    zH = (dp-c)/np.sqrt(2)
    zF = -(dp+c)/np.sqrt(2)
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def afc1_dprime_to_p(dp, c):
    zH = (dp-c)
    zF = -(dp+c)
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def dprime_to_p(dp, c):
    zH = -c+dp
    zF = -c
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def fake_dprime(h, f=None):
    if f is None:
        p = h
    else:
        p = (h + (1-f)) / 2.0
    z = ztop(p)
    dprime = z * np.sqrt(2)
    return dprime

