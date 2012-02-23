%for fitting a probit to data with an arbitrary lower asymptote
clear all; close all;format compact;clc
disp('******************************')
dataset=1;   %1 for p. 69 data,  2 for p. 94 data
if dataset==1,   % Data on p. 69 of Kingdom/Prins
    StimLevels=[.01 .03 .05 .07 .09 .11]   %see details on p. 69
    NumPos=[45 55 72 85 91 100]
    OutOfNum=ones(1,length(NumPos))*100 %number trials
    paramsInit=[.05 50]; %Thresh, slope & Lower Asymptote (latter is fixed for now)
elseif dataset==2,   %Data on p. 94 of Kingdom/Prins
    StimLevels=[-2 -1 0 1 2]
    NumPos= [2 3 3 3 4]
    OutOfNum=[4 4 4 4  4]
    paramsInit=[0 1];%Thresh, slope & Lower Asymptote (latter is fixed for now)
elseif dataset==3, %p 118 Example 1
    StimLevels=-2:1:2
    NumPos=[48 53 55 100 100]
    OutOfNum = 100.*ones(size(StimLevels))
    paramsInit=[0.5 4]
end
LowerAsymptote=.5;
ProbitOrLogit=2;%1,3 is probit, 2,4 is logit  (probit means cumulative normal)
        %for ProbitOrLogit<=2 z=(stim-p1)*p2   like Kingdom&Prins
        %for ProbitOrLogit>2 z=(stim-p1)/p2   like p2 being like JND in stim units
[params1,chisqLL]=fminsearch('ProbitLogit',paramsInit,[], StimLevels, NumPos, OutOfNum, ...
    LowerAsymptote, ProbitOrLogit,1);

%% now do the search
midstim=(StimLevels(end)+StimLevels(1))/2;
range=StimLevels(end)-StimLevels(1);
N=20;
inc=range*4/N
for iTh=1:N;
    for iSl=1:N;  %slope    horizontal
        Thresh(iTh)=midstim+iTh*inc-2*range;
        Slope(iSl)=2*(iSl-5)/range;
        paramsInit=[Thresh(iTh) Slope(iSl)];
        [params,chisqLL]=fminsearch('ProbitLogit',paramsInit,[], ...
            StimLevels, NumPos,OutOfNum, LowerAsymptote, ProbitOrLogit,1);
        error(iTh,iSl)=log((params(1)-params1(1))^2+(params(2)-params1(2))^2);
    end
end
contourf(Thresh,Slope,error');colorbar
xlabel('threshold');ylabel('slope')
title('log of the sum of square error of the 2 parameters')

