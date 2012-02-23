# for fitting a probit to data with an arbitrary lower asymptote
import numpy as np
import matplotlib.pyplot as plt
from pypsy.pf import ProbitLogit
import scipy.optimize as optimize
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

log.info('*'*30)

dataset = 1 #1 for p. 69 data,  2 for p. 94 data
numBootstraps = 10
plot_opt = 'both'
param_Init =  None
ProbitOrLogit=4;
#1 is probit, 2 is logit (probit means cumulative normal, prins slope) 
#3 is probit, 4 is logit (param2 becomes like JND, since it is reciprocated)

if len(sys.argv) > 1:
    dataset = int(sys.argv[1])
if len(sys.argv) > 2:
    numBootstraps = int(sys.argv[2])
if len(sys.argv) > 3:
    plot_opt = sys.argv[3]
if len(sys.argv) > 4:
    ProbitOrLogit = int(sys.argv[4])
if len(sys.argv) > 5:
    param_Init  = np.array([float(anarg) for anarg in sys.argv[5:]])

if dataset==1:   # Data on p. 69 of Kingdom/Prins
    StimLevels=np.array([.01, .03, .05, .07, .09, .11])   #see details on p. 69
    NumPos=np.array([45., 55, 72, 85, 91, 100])
    Ntrials=np.ones(6)*100; #number trials
    paramInit=np.array([.05, 50, .5, 0]); #JND, PSE, Lower Asymptote (latter is fixed for now)
elif dataset==2:   #Data on p. 94 of Kingdom/Prins
    StimLevels=np.array([-2, -1, 0, 1, 2]);
    NumPos= np.array([2., 3, 3, 3, 3.999])*1;
    NumPos= np.array([2., 3, 3, 3, 4])*1;
    Ntrials=np.array([4., 4, 4, 4,  4])*1;
    paramInit=np.array([0, 1, .5, 0]); #JND, PSE, Lower Asymptote (latter is fixed for now)
    #if 1==0:
    #    paramsFree=np.array([1, 1, 0, 0]);
    #    #PF=@PAL_Logistic;
    #    paramsValues, LL, exitflag, output=PAL_PFML_Fit(StimLevels, NumPos,
    #        Ntrials, paramInit, paramsFree, PF);
elif dataset==3:   # p118 of Kingdom/Prins / 4.1
    StimLevels=np.array([-2., -1, 0, 1, 2])   #see details on p. 69
    NumPos=np.array([48., 53, 55, 100, 100])
    Ntrials=np.ones(5)*100; #number trials
    paramInit=np.array([-1, 1, .5, 0]); #JND, PSE, Lower Asymptote (latter is fixed for now)

if param_Init is not None:
    paramInit = param_Init
pObs=NumPos/Ntrials;  #probability correct
params0=paramInit[0:2];
LowerAsymptote=paramInit[2];

log.info("running dataset #%i",dataset)
log.info ('*'*30)
log.info('from initial conditions')
LogLik, p0=ProbitLogit(params0, StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1)
log.info('LogLik of initial: %.4g' % (LogLik))
if plot_opt in ('both','pf'):
    plt.subplot(1,1,1)
    plt.plot(StimLevels, pObs,'*r')
    # Plot a smoother fitted function
    smoothrang = np.linspace(StimLevels[0], StimLevels[-1], 100 )
    LogLikX, smoothprob=ProbitLogit(params0, smoothrang, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,0)
    plt.plot(smoothrang, smoothprob, '--b', label='init conds');
    #plt.plot(StimLevels,p0,'.--b', label='initial conditions')
    plt.title('Fit data');
    plt.xlabel('Stimulus Intensity');
    plt.ylabel('Probability Correct')

## now do the search
#[params,chisqLL]=fminsearch('ProbitLogit',params0,[], StimLevels, NumPos, Ntrials,
#    LowerAsymptote, ProbitOrLogit,1)
def errfunc(*args):
    return ProbitLogit(*args)[0]

warn = 1
#while warn != 0:
out = optimize.fmin(errfunc, params0, args=(StimLevels, NumPos,
    Ntrials, LowerAsymptote, ProbitOrLogit, 1),full_output=1, retall=True,
    disp=0);
pfinal = out[0]  # Y
warn = out[4]; params0 = out[0]
pfinal = out[0]  # Y
searched_params = np.array( out[5] )

LogLikf, probExpect=ProbitLogit(pfinal, StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1)
if plot_opt in ('both','pf'):
    # Plot a smoother fitted function
    smoothrang = np.linspace(StimLevels[0], StimLevels[-1], 100 )
    LogLikX, smoothprob=ProbitLogit(pfinal, smoothrang, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,0)
    plt.plot(smoothrang, smoothprob, '-b', label='LL search');

    error=np.sqrt(probExpect*(1-probExpect)/Ntrials);
    plt.errorbar(StimLevels,probExpect,error, fmt=None, ecolor='b');


    plt.ylim(plt.ylim()[0],plt.ylim()[1]+.01)
    ##axis([-.1 16 0 1.05])
    #xlabel('Stimulus Intensity'); title('Fit based on likelihood search')
Nlevels=len(probExpect);
degfree=Nlevels-2.;    #predicted value of chisquare
ProbExact=1-special.gammainc(degfree/2., LogLikf/2.)
log.info('ProbExact = %.4g ' % ProbExact )
##[paramLSQ,chisqLSQ,fLSQ,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('ProbitLogit',
##    params,[],[],[], StimLevels, NumPos, Ntrials,LowerAsymptote, ProbitOrLogit,2);
if 1==1:   #for Log Likelihood
    if plot_opt in ('both','pf'):
        ax = plt.gca()
        kw =dict( ha='right',
                transform=ax.transAxes)
        results = 'chisqLL = %.2g\nJND = %.2g\nPSE = %.2g'
        plt.text(.95,.05 ,  results% (LogLikf,(1./pfinal[1]),pfinal[0]) , **kw)
        plt.legend(loc='best')
        plt.show()
    log.info('p[0] = %.4g ' % (pfinal[0]) )  #this give offset
    log.info('p[1] = %.4g ' % (pfinal[1]) ) #this prints out the inverse of slope
#    pass
#else:   #For chi square (not yet implimented
#    j=full(j);  #something about sparse matrices
#    cov=inv(j'*j);  #see Numerical Recipes for meaning of covariance matrix
#    SE=sqrt(diag(cov))';#Standard error
#    text(.2,.9,['JND = ' num2str(params(1))  ' +- ' num2str(SE(1)) ])
#    text(.2,.84,['PSE = ' num2str(params(2)) ' +- ' num2str(SE(2)) ])
#    print ['JND = ' num2str(params(1),3) ' +- ' num2str(SE(1),2)]
#    print ['PSE = ' num2str(params(2),3) ' +- ' num2str(SE(2),2)]
#end
log.info('best LogLik = %.4g' %(LogLikf))
#
if numBootstraps>1:
### Do parametric and nonparametric Monte Carlo simulations ('bootstraps')
    d = {}
    for iExpectOrObserved in [1,2]: #for parametric vs nonparametric
        if iExpectOrObserved==1:
            log.info('parametric bootstrap')
            prob=probExpect
        else:
            log.info('Nonparametric bootstrap')
            prob=pObs
        Nsim = numBootstraps;
        par = np.empty((Nsim,2))
        chisqLL2 = np.empty(Nsim)
        # XXX: we should look into using scipy.stats or even pymc to speed up
        # and leverage code others have written to do this for us
        # also - for our fitting routing - we can use something like:
        # scipy.stats.logistic.fit()
        for i in range(Nsim):    #MonteCarlo simulations to get standard errors of params
            N=Ntrials[0];
            NumPos=np.sum(np.random.rand(N,Nlevels)<np.ones((N,1))*prob, axis=0);#only for constant Ntrials
            #options = optimset('Display','off');
            #[par[i,:],chisqLL2(i)]=fminsearch('ProbitLogit',params,[], StimLevels,
            #    NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1);
            out = optimize.fmin(errfunc, pfinal, args=(StimLevels, NumPos,
                Ntrials, LowerAsymptote, ProbitOrLogit, 1),full_output=1, disp=0);
            sys.stdout.write('.')
            par[i,:]= out[0]
            chisqLL2[i] = out[1]
        sys.stdout.write('\n')
        SEParams=par.std(axis=0)
        covar=np.cov(par, rowvar=0)
        meanChi=chisqLL2.mean(axis=0)
        d[iExpectOrObserved] =  dict(
            meanParams=par.mean(axis=0)
            ,SEParams=SEParams
            ,covar=covar
            ,SqrtOfVar=np.sqrt(np.diag(covar)).T
            ,Correlation1=np.corrcoef(par, rowvar=0)
            ,Correlation2=covar[0,1]/np.prod(SEParams)
            ,meanChi=meanChi
            ,stdChi=chisqLL2.std(axis=0)
            ,SEchiPredicted=np.sqrt(2.*degfree)  #predicted SE of chisquareg
            ,pvalue_chisq=1.-special.gammainc(degfree/2., meanChi/2.)
            )
        if log.getEffectiveLevel() <= logging.INFO:
            dprint(d[iExpectOrObserved])

#end
### make contour plots

# XXX: it'd be nice to do contour plots *independent* of what dataset we're
# plotting (e.g. take the optimal solution we get, and then explore a
# parameter grid around that solution)
if plot_opt in ('both','contour'):
    xpts = 200
    ypts = 200
    if dataset==1:
        #TODO not sure what the grid ranges should be for DS 1
        p1 = np.linspace( -0.2, 1.0, xpts, endpoint=True)
        logp2=np.linspace( 0.1, 0.5, ypts, endpoint=True)
    elif dataset==2:
        p1 = np.linspace( -0.2, 1.0, xpts, endpoint=True)
        logp2=np.linspace( 0.1, 0.5, ypts, endpoint=True)
    elif dataset==3:
        p1 = np.linspace( -0.1, 1.0, xpts, endpoint=True)
        logp2=np.linspace( 0.0, 0.5, ypts, endpoint=True)

    LL=np.empty(p1.shape+logp2.shape)

    chimin = np.finfo(float).max
    for i1 in np.arange(len(p1)):
        for i2 in np.arange(len(logp2)):
            param=np.array([p1[i1],logp2[i2]]);
            [LogLik, p0]=ProbitLogit(param, StimLevels, NumPos,
                Ntrials, LowerAsymptote, ProbitOrLogit,1);
            LL[i1,i2]=LogLik;
            if np.isnan(LogLik)==False and LogLik<chimin:
                chimin = LogLik
    plt.figure(2)
    X=p1; # xxx: is this **0 supposed to be matrix exponentiation?
    Y=logp2;
    V=np.append( chimin,chimin+.125*2**np.arange(0,8.001)) # include chimin, then the others
    plt.contourf(X,Y,LL.T, V);
    #plt.contourf(X,Y,LL.T);
    plt.colorbar(format='%.4g')
    plt.xlabel('75% correct (a)')
    plt.ylabel('slope (b)')
    plt.show()
