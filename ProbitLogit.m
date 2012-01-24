function [LogLik, prob] = ProbitLogit(param,stim, Obs, N,  lower_asymptote, ProbitOrLogit, ChisqOrLL)
%stim are the stimulus levels
%Obs is the Observed data
%N is the number of trials at each level
%param(1) is the JND, param(2) is PSEs (in stim units)
upper_asymptote=1; 
z=(stim-param(1))*param(2);
if ProbitOrLogit==1,
    probFull=.5*(erf(z/sqrt(2))+1);%psychometric function going from 0 to 1
else,
    probFull=1./(1+exp(-z));
end
prob = lower_asymptote + (upper_asymptote-lower_asymptote)*probFull;
Expect=prob.*N;
if ChisqOrLL==1,
    LogLik = -2*sum((Obs.*log(Expect./(Obs+eps)) +(N-Obs).*log((N-Expect)./(N-Obs+eps))));
else Chisq= sum((Obs-Expect).^2./Expect./(1-prob));
end

%clg;hold off;
%plot(stim,Obs./N,'x',sim,prob,'-')
%xlabel('stimulus'); ylabel('prob correct')
