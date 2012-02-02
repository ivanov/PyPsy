import scipy.special # for erf
import numpy as np

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
#clg;hold off;
#plot(stim,Obs./N,'x',sim,prob,'-')
#xlabel('stimulus'); ylabel('prob correct')
