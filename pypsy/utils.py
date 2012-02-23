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
    """
    Convert z-score to percent

    Examples
    --------
    >>> ztop(0)
    0.5
    >>> ztop(1)
    0.84134474606854281
    >>> ztop(1) - ztop(-1)
    0.68268949213708563
    """
    p = 0.5 + 0.5*special.erf(z/np.sqrt(2))
    return p

def ptoz(p):
    """
    Percent to Z score

    Examples
    --------
    >>> ptoz(.5)
    0.0
    >>> ptoz(.841)
    0.99857627061565934
    >>> ptoz(.16)
    -0.99445788320975315
    """
    p = np.asanyarray(p)
    z = np.sqrt(2)*special.erfinv(p*2-1);
    return z

def afc2_p_to_dprime(ph, pf=None):
    """Using 2-AFC hit and false alarm rates, calculate the dprime, criterion,
    a logbias, as well as the percent correct

    Examples
    --------
    >>> r = afc2_p_to_dprime(.8, .1)
    >>> dprint(_dprime=r[0],criterion=r[1], lnB=r[2], percentCorrect=r[3])
    _dprime = 1.501
    criterion = 0.3111
    lnB = 0.467
    percentCorrect = 0.85


    Notes
    -----
    You can call also call this function with one parameter, in which case we
    stick to the same convention as Palamedes and use the first column for the
    hits, and the second column for the false alarms.
    """
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

def afc2_dprime_to_p(dprime, c):
    """Convert 2-AFC dprime and criterion c to hit and false alarm rates

    Examples
    --------
    >>> r = afc2_dprime_to_p(2, .1)
    >>> dprint(percent__hits=r[0], percent_false_alarms=r[1])
    percent__hits = 0.9104
    percent_false_alarms = 0.06878

    Notes
    -----
    Unlike the 1-AFC version, this function uses a sqrt(2) factor, see XXX for
    an explanation
    """
    zH = (dprime-c)/np.sqrt(2)
    zF = -(dprime+c)/np.sqrt(2)
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def afc1_dprime_to_p(dprime, c):
    """Convert 1-AFC dprime and criterion c to hit and false alarm rates

    Examples
    --------
    >>> r = afc1_dprime_to_p(2, .1)
    >>> dprint(percent__hits=r[0], percent_false_alarms=r[1])
    percent__hits = 0.9713
    percent_false_alarms = 0.01786
    """
    zH = (dprime-c)
    zF = -(dprime+c)
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def dprime_to_p(dprime, c):
    """Using a dprime and a criterion c, calculate the corresponding percent
    hit and percent false alarm rate.

    TODO: is this task / experimental paradigm specific? if it is, we should
    note what that is (1AFC? go-nogo?)

    Examples
    --------
    >>> r = dprime_to_p(2, .1)
    >>> dprint(percent__hits=r[0], percent_false_alarms=r[1])
    percent__hits = 0.9713
    percent_false_alarms = 0.4602
    """
    zH = -c+dprime
    zF = -c
    ph = ztop( zH)
    pf = ztop( zF)
    return ph,pf

def fake_dprime(h, f=None):
    """
    TODO: remove this?

    Examples
    --------
    """
    if f is None:
        p = h
    else:
        p = (h + (1-f)) / 2.0
    z = ztop(p)
    dprime = z * np.sqrt(2)
    return dprime

