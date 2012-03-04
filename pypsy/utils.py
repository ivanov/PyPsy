"some functions for PyPsy"
import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt
import pypsy.pf as pf
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def param_scatter(params, true_params, idxOfLapse=3, ax=None, needToConvert=False, plot_marginals=False):
    global axMargx, axMargy

    if needToConvert:
        fit_params_adj = np.array( [pf.convertToWichmann( p ) for p in params] )
    else:
        fit_params_adj = np.copy(params)

    true_pse = pf.fn_weibull_inv(0.5, true_params )
    true_slope = pf.fn_weibull_deriv( true_pse, true_params )

    laps_lowers = np.where( fit_params_adj[:,idxOfLapse] < 0.0001 )[0]
    laps_uppers = np.where( fit_params_adj[:,idxOfLapse] >= 0.0599)[0]
    laps_mids = np.where( np.all((fit_params_adj[:,idxOfLapse] >= 0.0001,fit_params_adj[:,idxOfLapse] < 0.0599),0) )[0]
    
    if ax==None:
        fig=plt.figure()
        ax = fig.add_subplot(1,1,1)

    if plot_marginals:
        divider = make_axes_locatable( ax)
        # create a new axes with a height of 1.2 inch above the axScatter
        axMargx = divider.new_vertical(1.2, pad=0.22, sharex=ax)
        # create a new axes with a width of 1.2 inch on the right side of the
        # axScatter
        # TODO: Figure out how to make the hists not overlap with the scatter
        axMargy = divider.new_horizontal(1.2, pad=0.22, sharey=ax)

    ax.plot( fit_params_adj[laps_lowers,0], fit_params_adj[laps_lowers,1], 'bx', label='$\lambda_{est}=0$' )
    ax.plot( fit_params_adj[laps_mids,0], fit_params_adj[laps_mids,1], 'g*', label='$\lambda_{est}\/{o.w.}$' )
    ax.plot( fit_params_adj[laps_uppers,0], fit_params_adj[laps_uppers,1], 'r.', label='$\lambda_{est}=0.6$' )
    ax.legend( loc='upper left')
    ax.axis('tight')

    x = [fit_params_adj[:,0].min(), fit_params_adj[:,0].max()]
    y = [fit_params_adj[:,1].min(), fit_params_adj[:,1].max()]

    ax.plot( [true_pse, true_pse], [fit_params_adj[:,1].min(), fit_params_adj[:,1].max()], 'k:' )
    ax.plot( [fit_params_adj[:,0].min(), fit_params_adj[:,0].max()], [true_slope, true_slope], 'k:' )
    ax.set_xlabel('pse (~%.4g)' % true_pse )
    ax.set_ylabel('slope@pse (~%.4g)' % true_slope)
    ax.loglog()

    mean0 = np.mean( fit_params_adj[laps_lowers], 0 )
    mean1 = np.mean( fit_params_adj[laps_mids], 0 )
    mean2 = np.mean( fit_params_adj[laps_uppers], 0 )

    xbins = 30
    ybins = 30
    if plot_marginals:
        fig.add_axes(axMargx)
        axMargx.hist( fit_params_adj[:,0], bins=xbins )
        fig.add_axes(axMargy)
        ax.axis('tight')
        axMargy.set_yscale("log", nonposy='clip')
        bins = 10**np.linspace( np.log10(y[0]), np.log10(y[1]), num=ybins )
        axMargy.hist( fit_params_adj[:,1],orientation='horizontal', bins=bins )
        oldticks = axMargy.get_xticklabels()
        [tick.set_rotation(-90) for tick in oldticks]
        ax.axis('tight')

    covarp=np.cov(fit_params_adj.T);
    sep=np.sqrt(np.diag(covarp))
    corp=covarp/(sep.T*sep);
    corrp=[corp[1,0],corp[idxOfLapse,0],corp[idxOfLapse,1]]

    covar0 = np.cov(fit_params_adj[laps_lowers].T); se0 = np.sqrt( np.diag( covar0))
    covar1 = np.cov(fit_params_adj[laps_mids].T); se1 = np.sqrt( np.diag( covar1))
    covar2 = np.cov(fit_params_adj[laps_uppers].T); se2 = np.sqrt( np.diag( covar2))

    tiny=8
    covartop = 0.3
    hite = 0.05
    ax.text( 0.05, covartop, r'$corr_{\alpha,\beta}$=%.4g' % corrp[0],
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.text( 0.05, covartop-hite, r'$corr_{\alpha,\lambda}$=%.4g' % corrp[1],
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
    ax.text( 0.05, covartop-2*hite, r'$corr_{\beta,\lambda}$=%.4g' % corrp[2],
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

    parmstop = 0.12
    hite2 = 0.04
    ax.text( 0.05, parmstop, '$\mu_{0.0}$=%s' % (mean_se_string(mean0, se0)),
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, size=tiny)
    ax.text( 0.05, parmstop-hite2, '$\mu$=%s' % (mean_se_string(mean1, se1)),
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, size=tiny)
    ax.text( 0.05, parmstop-hite2*2, '$\mu_{0.06}$=%s' % (mean_se_string(mean2, se2)),
            horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, size=tiny)

    return covarp, sep, corp, corrp

def mean_se_string(means,ses):
    thestr = ""
    for n in np.arange(len(means)-1):
        if n > 0:
            thestr += ', '
        thestr += '%.4g$\pm$%.4g' % (means[n],ses[n])
    return thestr


