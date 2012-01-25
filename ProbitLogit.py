import scipy.special # for erf
import numpy as np

def ProbitLogit(param,stim, Obs, N,  lower_asymptote, ProbitOrLogit, ChisqOrLL):
#stim are the stimulus levels
#Obs is the Observed data
#N is the number of trials at each level
#param(1) is the JND, param(2) is PSEs (in stim units)
# NOTE: above comment is wrong. first param is PSE, second is JND. drc + pi
    upper_asymptote=1; 
    z=(stim-param[0])*param[1];
    if ProbitOrLogit==1:
        probFull=.5*(scipy.special.erf(z/np.sqrt(2))+1);#psychometric function going from 0 to 1
    else:
        probFull=1/(1+np.exp(-z));
    prob = lower_asymptote + (upper_asymptote-lower_asymptote)*probFull;
    Expect=prob*N;
    if ChisqOrLL==1:
        LogLik = -2*sum((Obs*np.log(Expect/(Obs+np.finfo(float).eps))
            +(N-Obs)*np.log((N-Expect)/(N-Obs+np.finfo(float).eps))));
    else:
        Chisq= sum((Obs-Expect)**2/Expect/(1-prob));
    
    return [LogLik, prob] 
#clg;hold off;
#plot(stim,Obs./N,'x',sim,prob,'-')
#xlabel('stimulus'); ylabel('prob correct')
