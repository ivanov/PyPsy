# for fitting a probit to data with an arbitrary lower asymptote
import numpy as np
import matplotlib.pyplot as plt
from ProbitLogit import ProbitLogit
import scipy.optimize as optimize
import sys
import scipy.special as special

from utils import dprint

print('******************************')
dataset=1;   #1 for p. 69 data,  2 for p. 94 data
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
    if 1==0:
        paramsFree=np.array([1, 1, 0, 0]);
        #PF=@PAL_Logistic;
        paramsValues, LL, exitflag, output=PAL_PFML_Fit(StimLevels, NumPos,
            Ntrials, paramInit, paramsFree, PF);
pObs=NumPos/Ntrials;  #probability correct
params0=paramInit[0:2];
LowerAsymptote=paramInit[2];
ProbitOrLogit=2;#1 is probit, 2 is logit  (probit means cumulative normal)
disp('from initial conditions')
LogLik, p0=ProbitLogit(params0, StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1)
subplot(1,1,1)
plot(StimLevels, pObs,'*r')
plot(StimLevels,p0,'--b', label='initial conditions')
title('Fit data');
xlabel('Stimulus Intensity');ylabel('Probability Correct')

## now do the search
#[params,chisqLL]=fminsearch('ProbitLogit',params0,[], StimLevels, NumPos, Ntrials,
#    LowerAsymptote, ProbitOrLogit,1)
def errfunc(*args):
    return ProbitLogit(*args)[0]

warn = 1
while warn != 0:
    out = optimize.fmin(errfunc, params0, args=(StimLevels, NumPos,
        Ntrials, LowerAsymptote, ProbitOrLogit, 1),full_output=1); 
    pfinal = out[0]  # Y 
    warn = out[4]; params0 = out[0]
pfinal = out[0]  # Y

LogLikf, probExpect=ProbitLogit(pfinal, StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1)
plt.plot(StimLevels,probExpect,'-b', label='likelihood search')
error=np.sqrt(probExpect*(1-probExpect)/Ntrials);
plt.errorbar(StimLevels,probExpect,error,fmt=None, ecolor='b');
plt.ylim(plt.ylim()[0],plt.ylim()[1]+.01)
##axis([-.1 16 0 1.05])
##plt.text(.12,.55,'chisqLL = ' num2str(chisqLL,2)])
#xlabel('Stimulus Intensity'); title('Fit based on likelihood search')
Nlevels=len(probExpect);
degfree=Nlevels-2;    #predicted value of chisquare
ProbExact=1-special.gammainc(LogLikf/2,degfree/2) 
print('ProbExact = %.4g ' % ProbExact )
##[paramLSQ,chisqLSQ,fLSQ,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('ProbitLogit',
##    params,[],[],[], StimLevels, NumPos, Ntrials,LowerAsymptote, ProbitOrLogit,2);
if 1==1:   #for Log Likelihood
    text(.12,.5 , 'JND = %.2g' % (1./pfinal[1]) )
    text(.12,.45, 'PSE = %.2g' % pfinal[0])
    print('JND = %.4g ' % (1./pfinal[1]))  #this prints out the inverse of slope
    print('PSE = %.4g ' % (pfinal[0]) )  #this give offset
#    pass
#else:   #For chi square (not yet implimented
#    j=full(j);  #something about sparse matrices
#    cov=inv(j'*j);  #see Numerical Recipes for meaning of covariance matrix
#    SE=sqrt(diag(cov))';#Standard error
#    text(.2,.9,['JND = ' num2str(params(1))  ' +- ' num2str(SE(1)) ])
#    text(.2,.84,['PSE = ' num2str(params(2)) ' +- ' num2str(SE(2)) ])
#    disp(['JND = ' num2str(params(1),3) ' +- ' num2str(SE(1),2)])
#    disp(['PSE = ' num2str(params(2),3) ' +- ' num2str(SE(2),2)])
#end
print('chi square = %.2g' %(LogLikf))
plt.legend(loc='lower right')
plt.show()
#
### Do parametric and nonparametric Monte Carlo simulations ('bootstraps')
d = {}
for iExpectOrObserved in [1,2]: #for parametric vs nonparametric
    if iExpectOrObserved==1:
        disp('parametric bootstrap')
        prob=probExpect
    else:
        disp('Nonparametric bootstrap')
        prob=pObs
    Nsim=40;
    par = np.empty((Nsim,2))
    chisqLL2 = np.empty(Nsim)
    for i in range(Nsim):    #MonteCarlo simulations to get standard errors of params
        N=Ntrials[0];
        NumPos=sum(rand(N,Nlevels)<np.ones((N,1))*prob, axis=0);#only for constant Ntrials
        #options = optimset('Display','off');
        #[par[i,:],chisqLL2(i)]=fminsearch('ProbitLogit',params,[], StimLevels,
        #    NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1);
        out = optimize.fmin(errfunc, pfinal, args=(StimLevels, NumPos,
            Ntrials, LowerAsymptote, ProbitOrLogit, 1),full_output=1, disp=0); 
        sys.stdout.write('.')
        par[i,:]= out[0]
        chisqLL2[i] = out[1]
    sys.stdout.write('\n')
    SEParams=std(par, axis=0)
    covar=cov(par, rowvar=0)
    meanChi=mean(chisqLL2, axis=0)
    d[iExpectOrObserved] =  dict(
        meanParams=mean(par, axis=0)
        ,SEParams=SEParams
        ,covar=covar
        ,SqrtOfVar=sqrt(diag(covar)).T
        ,Correlation1=corrcoef(par, rowvar=0)
        ,Correlation2=covar[0,1]/prod(SEParams)
        ,meanChi=meanChi
        ,stdChi=std(chisqLL2, axis=0)
        ,SEchiPredicted=sqrt(2*degfree)  #predicted SE of chisquareg
        ,pvalue_chisq=1-special.gammainc(meanChi/2,degfree/2) 
        )
    dprint(d[iExpectOrObserved])
    #ProbExact=Stats2Prob('chi',meanChi, degfree, 0)
    #pvalue_chisq=1-special.gammainc(meanChi/2,degfree/2) 
    #print('pvalue_chisq = %.4g ' % pvalue_chisq )



#end
### make contour plots

# XXX: it'd be nice to do contour plots *independent* of what dataset we're
# plotting (e.g. take the optimal solution we get, and then explore a
# parameter grid around that solution)
if dataset==2:
    NumPos=[2, 3, 3, 3, 4];
    p1=np.arange(-2,2.001,.1);   #should be centered at params(1), extent given by SE(1)
    logp2=np.arange(-1.,1.,.1);
    LL=np.empty(p1.shape+logp2.shape)
    for i1 in range(len(p1)):
        for i2 in range(len(logp2)):
            param=np.array([p1[i1],10**logp2[i2]]);
            [LogLik, p0]=ProbitLogit(param, StimLevels, NumPos,
                Ntrials, LowerAsymptote, ProbitOrLogit,1);
            LL[i1,i2]=LogLik;
    chimin= LL.min()
    figure(2)
    X=p1; # xxx: is this **0 supposed to be matrix exponentiation?
    Y=logp2;
    V=chimin+.125*2**np.arange(0,8.001);
    plt.contourf(X,Y,LL.T,V);
    plt.colorbar()
    xlabel('75% correct (a)')
    ylabel('log slope (b)')

plt.show()
