"some functions for PyPsy"
import numpy as np
import scipy.special as special
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

