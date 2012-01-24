%for fitting a probit to data with an arbitrary lower asymptote
clear all; close all
disp('******************************')
dataset=1;   %1 for p. 69 data,  2 for p. 94 data
if dataset==1,   % Data on p. 69 of Kingdom/Prins
    StimLevels=[.01 .03 .05 .07 .09 .11]   %see details on p. 69
    NumPos=[45 55 72 85 91 100]
    Ntrials=ones(1,6)*100; %number trials
    paramInit=[.05 50 .5 0]; %JND, PSE, Lower Asymptote (latter is fixed for now)
elseif dataset==2,   %Data on p. 94 of Kingdom/Prins
    StimLevels=[-2 -1 0 1 2];
    NumPos= [2 3 3 3 3.999]*1;
    NumPos= [2 3 3 3 4]*1;
    Ntrials=[4 4 4 4  4]*1;
    paramInit=[0 1 .5 0]; %JND, PSE, Lower Asymptote (latter is fixed for now)
    if 1==0;
        paramsFree=[1 1 0 0];
        PF=@PAL_Logistic;
        [paramsValues LL exitflag output]=PAL_PFML_Fit(StimLevels, NumPos, ...
            Ntrials, paramInit, paramsFree, PF);
        return
    end
end
pObs=NumPos./Ntrials;  %probability correct
params0=paramInit(1:2);
LowerAsymptote=paramInit(3);
ProbitOrLogit=2;%1 is probit, 2 is logit  (probit means cumulative normal)
disp('from initial conditions')
[LogLik, p0]=ProbitLogit(params0, StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1)
subplot(1,2,1)
plot(StimLevels, pObs,'*r', StimLevels,p0,'--b')
title('Fit data based in initial guesses');
xlabel('Stimulus Intensity');ylabel('Probability Correct')

%% now do the search
[params,chisqLL]=fminsearch('ProbitLogit',params0,[], StimLevels, NumPos, Ntrials, ...
    LowerAsymptote, ProbitOrLogit,1)
stimPlot=[-1:.1:8];
[chisqLL, probExpect]=ProbitLogit(params,StimLevels, NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1) %for making plot
%[deviation, prob]=ProbitLogit(params, 1, stim, probs,lowerAsymp,optim,ProbLog) %to get prob
subplot(1,2,2);
plot(StimLevels,probExpect,'b',StimLevels,pObs,'r*');hold on
error=sqrt(probExpect.*(1-probExpect)./Ntrials);
errorbar(StimLevels,probExpect,error,'.');
%axis([-.1 16 0 1.05])
text(.12,.55,['chisqLL = ' num2str(chisqLL,2)])
xlabel('Stimulus Intensity'); title('Fit based on likelihood search')
Nlevels=length(probExpect);
degfree=Nlevels-2;    %predicted value of chisquare
ProbExact=Stats2Prob('chi',chisqLL, degfree, 0)
%[paramLSQ,chisqLSQ,fLSQ,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('ProbitLogit',...
%    params,[],[],[], StimLevels, NumPos, Ntrials,LowerAsymptote, ProbitOrLogit,2);
if 1==1,   %for Log Likelihood
    text(.12,.5 ,['JND = ' num2str(1/params(2),2) ])
    text(.12,.45,['PSE = ' num2str(params(1),2)])
    disp(['JND = ' num2str(1/params(2),4)])  %this prints out the inverse of slope
    disp(['PSE = ' num2str(params(1),4) ])  %this give offset
else   %For chi square (not yet implimented
    j=full(j);  %something about sparse matrices
    cov=inv(j'*j);  %see Numerical Recipes for meaning of covariance matrix
    SE=sqrt(diag(cov))';%Standard error
    text(.2,.9,['JND = ' num2str(params(1))  ' +- ' num2str(SE(1)) ])
    text(.2,.84,['PSE = ' num2str(params(2)) ' +- ' num2str(SE(2)) ])
    disp(['JND = ' num2str(params(1),3) ' +- ' num2str(SE(1),2)])
    disp(['PSE = ' num2str(params(2),3) ' +- ' num2str(SE(2),2)])
end
disp(['chi square = ' num2str(chisqLL,2)])

%% Do parametric and nonparametric Monte Carlo simulations ('bootstraps')
for iExpectOrObserved=1:2, %for parametric vs nonparametric
    if iExpectOrObserved==1,
        disp('parametric bootstrap')
        prob=probExpect
    else
        disp('Nonparametric bootstrap')
        prob=pObs
    end
    Nsim=400;
    for i=1:Nsim;    %MonteCarlo simulations to get standard errors of params
        N=Ntrials(1);
        NumPos=sum(rand(N,Nlevels)<ones(N,1)*prob);%only for constant Ntrials
        options = optimset('Display','off');
        [par(i,:),chisqLL2(i)]=fminsearch('ProbitLogit',params,[], StimLevels, ...
            NumPos, Ntrials, LowerAsymptote, ProbitOrLogit,1);
    end
    meanParams=mean(par)
    SEParams=std(par)
    covar=cov(par)
    SqrtOfVar=sqrt(diag(covar))'
    Correlation1=corrcoef(par)
    Correlation2=covar(1,2)/prod(SEParams)
    meanChi=mean(chisqLL2)
    stdChi=std(chisqLL2)
    SEchiPredicted=sqrt(2*degfree)  %predicted SE of chisquareg
    ProbExact=Stats2Prob('chi',meanChi, degfree, 0)
    
end
%% make contour plots
LL=[];
if dataset==2,
    NumPos=[2 3 3 3 4];
    p1=-2:.1:2;   %should be centered at params(1), extent given by SE(1)
    logp2=-1:.1:1;
    for i1= 1:length(p1),
        for i2=1:length(logp2);
            param=[p1(i1) 10.^logp2(i2)];
            [LogLik, p0]=ProbitLogit(param, StimLevels, NumPos, ...
                Ntrials, LowerAsymptote, ProbitOrLogit,1);
            LL(i1,i2)=LogLik;
        end
    end
    chimin=min(min(LL));
    figure(2)
    X=logp2'.^0*p1;
    Y=logp2'*p1.^0;
    V=chimin+.125*2.^[0:8];
    contourf(X,Y,LL',V);colorbar
    xlabel('75% correct (a)')
    ylabel('log slope (b)')
end