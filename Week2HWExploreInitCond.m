%for fitting a probit to data with an arbitrary lower asymptote
clear all; close all;format compact
disp('******************************')
dataset=1;   %1 for p. 69 data,  2 for p. 94 data
if dataset==1,   % Data on p. 69 of Kingdom/Prins
    StimLevels=[.01 .03 .05 .07 .09 .11]   %see details on p. 69
    NumPos=[45 55 72 85 91 100]
    Ntrials=ones(1,6)*100; %number trials
    params0=[.05 50]; %JND, PSE, Lower Asymptote (latter is fixed for now)
elseif dataset==2,   %Data on p. 94 of Kingdom/Prins
    StimLevels=[-2 -1 0 1 2];
    NumPos= [2 3 3 3 4]*1;
    Ntrials=[4 4 4 4  4]*1;
    params0=[0 1]; %JND, PSE, Lower Asymptote (latter is fixed for now)
end
LowerAsymptote=.5;
ProbitOrLogit=4;%1,3 is probit, 2,4 is logit  (probit means cumulative normal)
        %for ProbitOrLogit<3 z=(stim-p1)*p2   like Kingdom&Prins
        %for ProbitOrLogit>2 z=(stim-p1)/p2   like p2 being like JND in stim units
[params1,chisqLL]=fminsearch('ProbitLogit',params0,[], StimLevels, NumPos, Ntrials, ...
    LowerAsymptote, ProbitOrLogit,1);

%% now do the search
for iTh=1:40;
    for iSl=1:20;  %slope    horizontal
        if dataset==1,
            Thresh(iTh)=(iTh-5)/20;   %Threshold  vertical
            Slope(iSl)=(iSl-10)*10;
        else
            Thresh(iTh)=(iTh-20);   %Threshold  vertical
            Slope(iSl)=(iSl-10)*10;
        end
        params0=[Thresh(iTh) Slope(iSl)];
        [params,chisqLL]=fminsearch('ProbitLogit',params0,[], StimLevels, NumPos, Ntrials, ...
            LowerAsymptote, ProbitOrLogit,1);
        error(iTh,iSl)=log((params(1)-params1(1))^2+(params(2)-params1(2))^2);;
    end
end
contourf(Thresh,Slope,error');colorbar
xlabel('threshold');ylabel('slope')
title('log of the sum of square error of the 2 parameters')

