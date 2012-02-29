"""Psychometric function related code

Currently, only Probit and Logit are implemented, but the `fitpf` is generic
enough to take anything you might pass it as errfunc. If that doesn't work,
it's considered a bug, so please file it as such.
"""

import numpy as np
import scipy.special # for erf
import scipy.optimize as optimize
import scipy.misc as misc

def ProbitLogit(param,stim, Obs, N,  lower_asymptote, ProbitOrLogit, ChisqOrLL):
    """
    Given all of the parameters, return either (depending on the value of
    ChisqOrLL) the log-likelihood or the chi-squared value along with the
    probability for the data under the model

    Parameters
    ----------
    param: array of length 2
        param[0] is the PSE (in stim units), param[1] is the either the slope or
        the JND, depending on the value of ProbitOrLogit
    stim : array
        The stimulus levels
    Obs : array
        The Observed data: number of correct responses at each level
    N : int
        The number of trials at each level
    lower_asymptote : float
        The 'minimum' used when rescaling (e.g.: 0.5 for binary outcomes)
    ProbitOrLogit : int (1,2,3, or 4)
        1 is probit, 2 is logit (probit means cumulative normal, prins slope)
        3 is probit, 4 is logit (param[1] becomes like JND, since it is
        reciprocated)
    ChisqOrLL : int (1,2, or 0)
        1 is Log-Likelihood, 2 is Chi-Squared value

    Returns
    -------
    LogLik, prob : float,float 
        The log-likelihood or the chi-squared value along with the
        probability for the data under the model

    Notes
    -----
    Original Matlab code *wrongly* stated that:
         "param(1) is the JND, param(2) is PSEs (in stim units)"
    We have corrected that here - drc + pi
    """
    upper_asymptote=1.0; 

    tparam = param.copy()

    if ProbitOrLogit<1 or ProbitOrLogit > 4:
        raise ValueError("Don't know what to do when ProbitOrLogit = %i"%ProbitOrLogit)
    if ProbitOrLogit>2:
        tparam[1]=1./tparam[1]; #param(2) is now like JND rather than slope (in stim units)
        ProbitOrLogit=ProbitOrLogit-2; #To allow it to be like before

    z=(stim-tparam[0])*tparam[1];

    if ProbitOrLogit==1:
        probFull=.5*(scipy.special.erf(z/np.sqrt(2))+1);#psychometric function going from 0 to 1
    else:
        probFull=1/(1+np.exp(-z));
    prob = lower_asymptote + (upper_asymptote-lower_asymptote)*probFull;

    if ChisqOrLL>0:
        Expect=prob*N

    if ChisqOrLL==1:
        LogLik = -2*sum((Obs*np.log(Expect/(Obs+np.finfo(float).eps))
            +(N-Obs)*np.log((N-Expect)/(N-Obs+np.finfo(float).eps))));
    elif ChisqOrLL==2:
        LogLik= sum((Obs-Expect)**2./Expect/(1.-prob))
    elif ChisqOrLL==0:
        LogLik=0;  #for plotting
    else:
        raise ValueError("Only know ChisqOrLL==0,1,2, not %d"%ChisqOrLL)

    return LogLik, prob

def errfunc(*args):
    "simple wrapper of ProbitLogit for use with fmin"
    return ProbitLogit(*args)[0]

def fitpf(params0, StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,
        output_param_search=False, errfunc=errfunc):
    """Fit a psychometric function.

    The default function to minimize is ProbitLogit.

    XXX: make this more generic to fit other functions, so we don't have to
    pass inelegant "ProbitOrLogit"-type of parameters
    """
    warn = 1
    #while warn != 0:
    out = optimize.fmin(errfunc, params0, args=(StimLevels, NumPos,
        Ntrials, LowerAsymptote, ProbitOrLogit, 1),full_output=1, retall=True,
        disp=0, xtol=1e-6, ftol=1e-6,maxiter=400*len(params0),maxfun=400*len(params0));
    pout = out[0]  # Y
    warn = out[4]; params0 = out[0]
    pfinal = out[0]  # Y
    searched_params = np.array( out[5] )

    if output_param_search:
        return pout, searched_params
    else:
        return pout

class experiment():
    def __init__( self, levels, Ntrials, Ncorr=None, pcorr=None ):
        self.levels = np.array(levels)
        self.Ntrials = np.array(Ntrials, dtype='int')
        self.Ncorr = np.array(Ncorr)

        if pcorr is None:
            self.pcorr = Ncorr/np.array(Ntrials, dtype='float')
        else:
            self.pcorr = pcorr

            if Ncorr is None:
                self.simulate()

    def simulate(self):
        self.Ncorr = np.random.binomial( self.Ntrials, self.pcorr )
        return self.Ncorr

# For use inside pf_generic:
def fn_probit(x):
    return .5*(scipy.special.erf(x/np.sqrt(2.))+1.);
def fn_logit(x):
    return 1./(1.+np.exp(-x));

def fn_logistic(x, params):
    # params[0]==mu, params[1]==theta
    return 1./(1.+np.exp(-(x-params[0])/params[1]));

# TODO: Put these in a class abstracting a fn, params,
# deriv and inverse?

# TODO: Validation!!
# weibull([10,3,0,0])
#fw2 = [4.7231,7.0918,7.9939,9.7128,10.6386,13.2050]
# From Wichmann&Hill/Prins:
# pse = weibull_inv(0.5, [10,3,0,0])
# pse ~ 8.85
# weibull( pse, [10,3,0,0] ) ~ 0.5
# weibull_deriv( pse, [10,3,0,0] ) ~ 0.1175 (Wichmann says '0.118' in paper)
def fn_weibull(x, params):
    [alpha,beta,gamma,lambd]=params
    y = gamma + (1 - gamma - lambd)*(1 - np.exp(-(x/alpha)**beta))
    return y
def fn_weibull_inv(x, params):
    [alpha,beta,gamma,lambd]=params
    c = (x-gamma)/(1.0-gamma-lambd)
    temp1 = -np.log(1.0-c)
    # Really just want alpha*temp**(1/beta), but
    # Python doesn't like to take a fractional power
    # of a negative number (MATLAB seems okay with that)
    y = alpha*np.exp((1.0/beta)*np.log(temp1))
    return y
def fn_weibull_deriv( x, params):
    [alpha,beta,gamma,lambd]=params
    y = (1-gamma-lambd)*np.exp(-(x/alpha)**beta)*(x/alpha)**(beta-1)*beta/alpha
    return y

def errfunc_OO(*args):
    "simple wrapper for PF class fitting using fmin"
    "0: (params that fmin is exploring...)"
    "1: self"
    "2: data"
    "3: which stat to use (1=LL, 2=X2, 3=L_TreutweinStrasburger)"
    obj=args[1]

    for learnable_param_idx,target_param_idx in enumerate(obj.params_free):
        obj.params[target_param_idx] = args[0][learnable_param_idx]

    return obj.eval_gof(args[2])[args[3]]

class pf_generic():
    # TODO: make params a dict ?
    PARAM_PSE=0
    PARAM_SPREAD=1
    PARAM_LOWER=2
    PARAM_UPPER=3
    PARAM_INVERT_SLOPE=4
    def __init__(self, fn, params, params_free=np.array([0,1]), invert_slope=False):
        self.fn = fn
        self.params = params
        self.params_free = params_free
        self.invert_slope = invert_slope

    # TODO: memoize results of next two functions?
    def eval( self, x):
        if self.invert_slope:
            spread = 1.0/self.params[self.PARAM_SPREAD]
        else:
            spread = self.params[self.PARAM_SPREAD]
        pse = self.params[self.PARAM_PSE]
        # rescale ordinate
        #eval_x=(x-pse)*spread
        self.probs = self.fn( x, self.params )
        #self.probs = self.params[self.PARAM_LOWER] + (1.0 - self.params[self.PARAM_UPPER]-self.params[self.PARAM_LOWER])*probs;
        return self.probs

    def eval_gof( self, data):
        probs = self.eval(data.levels)
        expected = probs * data.Ntrials
        # Compute various goodness of fit statistics
        LL = -2*sum((data.Ncorr*np.log(expected/(data.Ncorr+np.finfo(float).eps))
            +(data.Ntrials-data.Ncorr)*np.log((data.Ntrials-expected)/(data.Ntrials-data.Ncorr+np.finfo(float).eps))));
        X2 = sum((data.Ncorr-expected)**2./expected/(1.-probs))

        # Adding eps to avoid log(0) Nans
        self.prinsNLL = -sum( data.Ncorr*np.log(probs+np.finfo(float).eps)+(data.Ntrials-data.Ncorr)*np.log(1.-probs+np.finfo(float).eps) )

        # Treutwein/Strasburger 1999 Eq 6 (likelihood of the data)
        L_ts = 2**(sum( data.Ntrials ))
        LL_ts = 0.0
        #L_ts = 1.0
        for level in np.arange( len(data.levels) ):
            # TODO: is right to use observed data or function values?: next two lines can chg to try fitted
            thisN = data.Ntrials[level]
            thisCorr = data.Ncorr[level]
            L_ts *= misc.comb( thisN, thisCorr ) * (probs[level]**thisCorr) * (1.0 - probs[level])**(thisN-thisCorr) 
            LL_ts += np.log(misc.comb( thisN, thisCorr )) + thisCorr*np.log(probs[level]) +np.log(1.0 - probs[level])*(thisN-thisCorr) 

        #TODO: This is how Prins' clamps the lapse.  Parameterize.
        if (self.params[self.PARAM_UPPER] < 0) or (self.params[self.PARAM_UPPER] > 0.05):
            self.prinsNLL=np.inf

        return probs,LL,X2,L_ts,LL_ts,self.prinsNLL

    def fitpf(self, params0, data, output_param_search=False, errfunc=errfunc_OO, which_stat_to_min=5):
        """Fit a psychometric function.
        """
        warn = 1
        #while warn != 0:
        if params0[2]==0.0:
            params0[2]=0.00025 # like in PAL_minimize
        print params0
        out = optimize.fmin(errfunc_OO, params0, args=(self,data,which_stat_to_min), full_output=1, retall=True, disp=0);
        pout = out[0]  # Y
        warn = out[4]; params0 = out[0]
        pfinal = out[0]  # Y
        self.searched_params = np.array( out[5] )
    
        if output_param_search:
            return pout, self.searched_params
        else:
            return pout

    # This may not belong here. Currently unused:
    def findMaxGrid(self, x, param_grid):
        global values

        if getattr( x, '__iter__', False):
            pass
        else:
            x = [x]

        # OTHER method:
        #try:
            #iterator = iter(x)
        #except:
            #x = [x]

        values = np.zeros( np.concatenate( ([len(x)], [len(p) for p in param_grid] ) ) )
        for n0,p0 in enumerate(param_grid[0]):
            for n1,p1 in enumerate(param_grid[1]):
                for n2,p2 in enumerate(param_grid[2]):
                    values[:,n0,n1,n2] = np.log(self.fn( x, [p0,p1,0.5,p2]))

        return values

class pf_stan(pf_generic):
    # Only things different: The eval can invert the slope, and uses Stan's param
    # interpretation method, which is nonstandard. (vs. Prins/Wichmann)
    # In practice this can probably be removed and moved into
    # one of the pf_functions, I think...
    # This only works with logit and (maybe) probit.
    def __init__(self, fn, params, params_free=np.array([0,1]), invert_slope=True):
        pf_generic.__init__(self,fn,params,params_free,invert_slope)

    def eval( self, x):
        if self.invert_slope:
            spread = 1.0/self.params[self.PARAM_SPREAD]
        else:
            spread = self.params[self.PARAM_SPREAD]
        pse = self.params[self.PARAM_PSE]
        # rescale ordinate
        eval_x=(x-pse)*spread
        probs = self.fn( eval_x ) 
        self.probs = self.params[self.PARAM_LOWER] + (self.params[self.PARAM_UPPER]-self.params[self.PARAM_LOWER])*probs;
        return self.probs

