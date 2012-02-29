# for fitting a probit to data with an arbitrary lower asymptote
import numpy as np
import matplotlib.pyplot as plt
from pypsy.pf import ProbitLogit, fitpf
import pypsy.pf as pf
import sys
import scipy.special as special

from pypsy.utils import dprint

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

numtrials = 120
numlevels = 6
alpha = 10.0
beta = 3.0
lapse = 0.02
levels_frac = np.array([.1,.3,.4,.6,.7,.9]) # Prins' s2
params=[alpha,beta,0.5,lapse]

# Model:
fit1 =  pf.pf_generic( pf.fn_weibull, params, [0,1])

# Data:
levels = pf.fn_weibull_inv( levels_frac, [alpha,beta,0,0] )
trials_arr = np.tile( numtrials/numlevels, numlevels)
data=pf.experiment(levels, trials_arr, pcorr=pf.fn_weibull( levels, [alpha,beta,0.5,0]) )
data.simulate()

params0 = [9.,2.]

fit1.fitpf( params0, data )
gof = fit1.eval_gof( data )

log.info( str(gof) )
#log.info( "probs same? %s" % str( np.all( gof[0]==probExpect) ) )
#log.info( "LL same? %s" % str(gof[1]==LogLikf) )
#log.info( "X2 same? %s" % str(gof[2]==LogLikX2) )
log.info( "fitted params: %s " % str(testpf.params ) )
