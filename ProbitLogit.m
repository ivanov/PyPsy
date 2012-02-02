function [LogLik, prob] = ProbitLogit(param,stim, Obs, N, ...
    lower_asymptote,ProbitOrLogit, ChisqOrLL)  %fitting Probit or Logit
%Outputs:  prob, the Expected probability is the Psychometric Function (PL)
%LogLik is log likelihood in with chi square normalization.
%stim are the stimulus levels
%Obs is the Observed data (number correct)
%N is a vector with the number of trials at each level
%param(1)= the horizontal position of the PL (sometimes called PSE)
%param(2)=PL slope(if ChisqOrLL<2) otherwise= 1/slope in stimulus units. 
%lower_asymptote is a fixed lower asymptote
%ProbitOrLogit=1 is probit, =2 is logit, 
%ProbitOrLogit=3&4 is like 1&2 but use JND =1/slope to indicate steepness
%ChisqOrLL =1 likelihood, 2=chisq,  0 for plotting
upper_asymptote=1; 
if ProbitOrLogit>2, 
    param(2)=1/param(2); %param(2) is now JND (in stim units), not slope
    ProbitOrLogit=ProbitOrLogit-2;%To allow it to be like before
end
z=(stim-param(1))*param(2);  %(but look at line 12 for possible inverse)
if ProbitOrLogit==1,    %This is probit formula
    probFull=.5*(erf(z/sqrt(2))+1);%psychometric function going from 0 to 1
else,                   %This is logit formula
    probFull=1./(1+exp(-z));
end
prob = lower_asymptote + (upper_asymptote-lower_asymptote)*probFull;
if ChisqOrLL>0, Expect=prob.*N; end
if ChisqOrLL==1,
    LogLik = -2*sum((Obs.*log((Expect+eps)./(Obs+eps)) + ...
        (N-Obs).*log((N-Expect)./(N-Obs+eps))));
elseif ChisqOrLL==2,LogLik= sum((Obs-Expect).^2./Expect./(1-prob+eps));
elseif ChisqOrLL==0, LogLik=0;  %for plotting
end
%if isnan(LogLik), [Expect-Obs Obs], end
%clg;hold off;
%plot(stim,Obs./N,'x',sim,prob,'-')
%xlabel('stimulus'); ylabel('prob correct')
