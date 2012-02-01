%for fitting a probit to data with an arbitrary lower asymptote
clear all; close all
disp('******************************')
dataset=2;   %1 for p. 69 data,  2 for p. 94 data
if dataset==1,   % Data on p. 69 of Kingdom/Prins
    StimLevels=[.01 .03 .05 .07 .09 .11]   %see details on p. 69
    NumPos=[45 55 72 85 91 100]
    OutOfNum=ones(1,length(NumPos))*100 %number trials
    paramsInit=[.05 50 .5 0]; %JND, PSE, Lower Asymptote (latter is fixed for now)
elseif dataset==2,   %Data on p. 94 of Kingdom/Prins
    StimLevels=[-2 -1 0 1 2]
    NumPos= [2 3 3 3 4]
    OutOfNum=[4 4 4 4  4]
    paramsInit=[0 1 .5 0]; %JND, PSE, Lower Asymptote (latter is fixed for now)
elseif dataset==3,
    StimLevels=-2:1:2;
    NumPos=[48 53 55 100 100];
    OutOfNum = 100.*ones(size(StimLevels));
    paramsInit=[0.5 4 0.5 1]; %Initial, visually derived guesses for psyFunc
end
pObs=NumPos./OutOfNum;  %probability correct
params0=paramsInit(1:2);
LowerAsymptote=paramsInit(3);
ProbitOrLogit=2;%1 is probit, 2 is logit  (probit means cumulative normal)
disp('Get fit from initial conditions')
inc=(StimLevels(end)-StimLevels(1))/100;
StimPlot= StimLevels(1):inc:StimLevels(end);
[LogLik, p0]=ProbitLogit(params0, StimPlot, [], [], LowerAsymptote, ProbitOrLogit,0);
subplot(1,2,1)
plot(StimLevels, pObs,'*r', StimPlot,p0,'--b')
title('Fit data based in initial guesses');
xlabel('Stimulus Intensity');ylabel('Probability Correct')
%% Do the Palamedes program
disp('****The Palamedes program*******')
paramsFree=[1 1 0 0];
PF=@PAL_Logistic;
[paramsValues LL exitflag output]=PAL_PFML_Fit(StimLevels, NumPos, ...
    OutOfNum, paramsInit, paramsFree, PF)
%% Now do the search
disp('****Klein"s program*****')
ChisqOrLL=1; % 1-maximize likelihood like Prins,  2-minimize chi square

[params,chisqLL]=fminsearch('ProbitLogit',params0,[], StimLevels, NumPos, OutOfNum, ...
    LowerAsymptote, ProbitOrLogit,ChisqOrLL)
[chisqLL, probExpect]=ProbitLogit(params,StimLevels, NumPos, OutOfNum, LowerAsymptote, ProbitOrLogit,1); %for making plot
[chisqLL, probPlot]=ProbitLogit(params,StimPlot, [], [], LowerAsymptote, ProbitOrLogit,0); %for making plot
subplot(1,2,2);
plot(StimLevels,probExpect,'b.',StimLevels,pObs,'r*',StimPlot,probPlot,'b');hold on
error=sqrt(probExpect.*(1-probExpect)./OutOfNum);
errorbar(StimLevels,probExpect,error,'.');
text(StimLevels(end)*.9,.55,['chisqLL = ' num2str(chisqLL,2)])
xlabel('Stimulus Intensity'); title('Fit based on likelihood search')
Nlevels=length(probExpect);
degfree=Nlevels-2;    %predicted value of chisquare
pValueForChisq=1-gammainc(chisqLL/2,degfree/2)
if 1==1,   %If fminsearch was used so SE not calculated
    text(StimLevels(end)*.9,.5 ,['JND = ' num2str(1/params(2),2) ]) %for the plot
    text(StimLevels(end)*.9,.45,['PSE = ' num2str(params(1),2)])
    disp(['JND = ' num2str(1/params(2),4) ' This is inverse of Palamedes slope'])
    disp(['PSE = ' num2str(params(1),4) ' This is offset for 1AFC' ])
else   %To do nonlinear regression using fminsearch
    [paramLSQ,chisqLSQ,fLSQ,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('ProbitLogit',...
        params,[],[],[], StimLevels, NumPos, OutOfNum,LowerAsymptote, ProbitOrLogit,2);
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
for iExpectOrObserved=1:3, %for parametric vs nonparametric
    Nsim=400;   %this is the Prins B parameter for number of bootstraps
    if iExpectOrObserved==3,  %for Pamamedes version
        disp('***parametric bootstrap Palamedes***')
        [SD paramsSim LLSim converged]=PAL_PFML_BootstrapParametric( ...
            StimLevels, OutOfNum, [params .5 0], paramsFree,Nsim, PF) ;
        SD_Pal1=SD
        SD_Pal2=std(paramsSim)
    else %for Klein version
        disp('***Right now this bootstrap only works for equal #trials at all levels***')
        if iExpectOrObserved==1,
            disp('***parametric bootstrap Klein***')
            prob=probExpect
        elseif iExpectOrObserved==2,
            disp('***Nonparametric bootstrap Klein***')
            prob=pObs
        end
        for i=1:Nsim;    %MonteCarlo simulations to get standard errors of params
            N=OutOfNum(1);
            NumPos=sum(rand(N,Nlevels)<ones(N,1)*prob);%only for constant OutOfNum
            options = optimset('Display','off');
            [par(i,:),chisqLL2(i)]=fminsearch('ProbitLogit',params,[], StimLevels, ...
                NumPos, OutOfNum, LowerAsymptote, ProbitOrLogit,1);
        end
        disp(['meanParams=' num2str(mean(par),3)])
        SEParams = std(par)
        covar=cov(par)
        disp(['SqrtOfVar(double check)=' num2str(sqrt(diag(covar))')])
        Correls=corrcoef(par);
        disp(['correlation from corrcoef=' num2str(Correls(1,2))])
        disp(['correlation from covariance (double check)=' num2str(covar(1,2)/prod(SEParams))])
        meanChi=mean(chisqLL2)
        disp(['stdChi=' num2str(std(chisqLL2),3)])
        disp(['SEchiPredicted=' num2str(sqrt(2*degfree),3)])  %predicted SE of chisquareg
        disp(['pvalue_chisq=' num2str(1-gammainc(meanChi/2,degfree/2),3)])  %p-value for chisq
    end
end
%% make contour plots
LL=[];
if dataset==1, NumPos=[45 55 72 85 91 100]; %needed after the bootstrap randomizing
elseif dataset==2, NumPos=[2 3 3 3 4];
elseif dataset==3, [48 53 55 100 100];
end
p1=-2:.1:2;   %should be centered at params(1), extent given by SE(1)
logp2=-1:.1:2;
for i1= 1:length(p1),
    for i2=1:length(logp2);
        param=[p1(i1) 10.^logp2(i2)];
        [LogLik, p0]=ProbitLogit(param, StimLevels, NumPos, ...
            OutOfNum, LowerAsymptote, ProbitOrLogit,1);
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