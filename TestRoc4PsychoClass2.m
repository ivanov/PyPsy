clear all;format compact;clc; close all
disp('****************Rating STD analysis*******************')
LLorChisq=2;  %=2 for chi square minimization that gives SE and correl.
datatype=1; 
if datatype==1,   %4x4 data with d'=1,2,3
    data=[84 50 16  3;13 34 34 13;3 13 34 34;0  3 16 50]
 %Note the 84 is #trials for blank that are BELOW the 1st criterion
    params0= [-1 0 1 1 2 3]*2;%3 criteria and 3 dprimes (try other initial guesses)
elseif datatype==2,%2x4 data with 4 responses
    data=[84 50 16  3;16 50 84 97]
    params0= [0  1 2 3]; %1 criterion and 3 dprimes
end
info(1)=0;          %info (1) and (2) are old capabilities of SetRocSimpler 
info(2)=0;          %I forget what info(1) and (2) are for now
info(3)=LLorChisq;  %LLorChisq
info(4)=1;          %dprimeType  =1 for individual d' being free
info(5)=0;          %offset for some d' options
ROC4PsychoClass(data, params0, info)
