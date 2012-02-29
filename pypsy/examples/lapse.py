# This file closely follows Prins' 'Figure3Scrutinize.m', which is an attempt
# to reproduce Wichmann & Hill 2001 parameter search, especially with regard to
# lapse rate.
import numpy as np
import matplotlib.pyplot as plt
import pypsy.pf as pf
import time

import logging
logging.basicConfig()
log = logging.getLogger(__file__)

# change .INFO to .WARNING to supress output
log.setLevel(logging.INFO)
#log.setLevel(logging.DEBUG) # uncomment this get really verbose

try:
    # define this variable and set it to True in your namespace
    # before you do a `run -i hw1.py` or `run -i exploreHw1.ipy` to silence
    # the logging output
    _be_quiet
except NameError:
    pass
else:
    if _be_quiet == True:
        log.setLevel(logging.WARNING)

#fw2 = [4.7231,7.0918,7.9939,9.7128,10.6386,13.2050]
#s2 levels should ~= this, which is copied from Prins' MATLAB code

nsims = 30
doPlots = False
gridGrain = 150
numtrials = 120
numlevels = 6
trials_per = 120/numlevels
alpha = 10.0
beta = 3.0
lapse = 0.00
levels_frac = np.array([.1,.3,.4,.6,.7,.9]) # Prins' s2
params2=[alpha,beta,0,0] # Wichmann/Prins use this to gen some initial values
params=[alpha,beta,0.5,lapse]

alphas = np.logspace( *np.log10( pf.fn_weibull_inv(np.array([0.01,0.99]), params2) ), num=gridGrain )
betas = np.logspace( -1.0, 2.0, gridGrain)
gamma = 0.5
lambdas = np.arange( 0, .06+np.finfo(float).eps, 0.005)

# Model:
fit1 =  pf.pf_generic( pf.fn_weibull, params, [0,1,3])

# Data:
levels = pf.fn_weibull_inv( levels_frac, params2  )
trials_arr = np.tile( trials_per, numlevels)
data=pf.experiment(levels, trials_arr, pcorr=pf.fn_weibull( levels, [alpha,beta,0.5,0]) )
data.simulate()

# Rebuild the huge grid of initial LogPs if we need to
# Takes about 16.5 secs on my Thinkpad x60
#
# Like Prins, we keep it here for efficiency, but is arguable.
try:
    if logpcorr[0,0,0,0] != 0:
            pass
except NameError:
    tic=time.time()
    log.info( 'Building huge pcorr grid for search (happens once per session).')
    param_grid = [alphas, betas, lambdas]
    param_dims = [len(p) for p in param_grid]

    # To put levels first:
    #logpcorr = np.zeros( np.concatenate( ([len(levels)], [len(p) for p in param_grid] ) ) )
    #logpincorr = np.zeros( np.concatenate( ([len(levels)], [len(p) for p in param_grid] ) ) )

    # To put levels last (eases some matrix arith)
    logpcorr = np.zeros( np.append( param_dims, [len(levels)] ) )
    logpincorr = np.zeros( np.append( param_dims, [len(levels)] ) )
    # TODO: Speed me up!! Life is too short...
    for n0,p0 in enumerate(param_grid[0]):
        for n1,p1 in enumerate(param_grid[1]):
            for n2,p2 in enumerate(param_grid[2]):
                logpcorr[n0,n1,n2,:] = pf.fn_weibull( levels, [p0,p1,0.5,p2])

    logpincorr = np.log(1.0-logpcorr)
    # May have some log(0) to fix, they may cause us problems with min/max
    logpincorr[ np.isinf( logpincorr) ] = np.log(np.finfo(float).eps )
    logpcorr = np.log(logpcorr)
    toc=time.time()

# NumPos = PAL_PF_SimulateObserver, etc...
NumPos = np.array([11,12,13,16,18,20]) # to confirm is working--enter simulation vals from Prins
NumPos = np.array([11,12,15,17,17,19]) # to confirm is working--enter simulation vals from Prins
data = pf.experiment( levels, trials_arr, NumPos )

if doPlots:
    smoothrang = np.linspace( levels[0], levels[-1], 30 )
    plt.ion()
    plt.figure()
    lplot=plt.subplot(1,2,1)
    rplot=plt.subplot(1,2,2)

fit_params = np.zeros( (nsims,len(params)) )
fit_LLs = np.zeros( nsims )
fit_data = np.zeros( (nsims, len(levels)) )
for asim in np.arange(nsims):
    if (asim%10)==0:
        log.warning( '%d/%d' % (asim,nsims))
    NumPos = data.simulate()
    LLspace = np.inner(logpcorr,NumPos) + np.inner(logpincorr, trials_arr-NumPos)
    maxidx = np.unravel_index( LLspace.argmax(), param_dims)
    p0 = [ alphas[maxidx[0]], betas[maxidx[1]], lambdas[maxidx[2]] ]
    fit1.fitpf( p0, data )
    gof = fit1.eval_gof( data )
    fit_params[asim] = fit1.params
    fit_LLs[asim] = fit1.prinsNLL
    fit_data[asim] = NumPos
    log.info( str(gof) )
    log.info( "fitted params: %s " % str(fit1.params ) )
    if doPlots:
        lplot.clear()
        lplot.contour( alphas, betas, np.max(LLspace,2))
        rplot.clear()
        rplot.plot( levels, NumPos, 'o' )
        rplot.plot( smoothrang, trials_per*fit1.eval( smoothrang), 'k-',lw=3 )
        plt.ylim( 5,22)
        plt.show()
        plt.draw()
