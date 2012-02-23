%Signal Detection Theory   Chapter 6
clear all;close all;  format compact; clc
dataset=1;
if dataset==1, %Data from p. 163 of Kingdom/Prins
    pHpF=[.2 .2; .6 .2; .7 .2; .8 .2; .84 .16]'
    x=[0 2:5];
elseif dataset==2, %from Exercise 2 p. 188
    pHpF=[.61 .69 .79 .88 .97 .99; .53 .42 .33 .18 .06 .03]
    x=[1:6];
end
pH=pHpF(1,:); pF=pHpF(2,:);     %hit and false alarm probabilities
%% Klein Method
disp('***      Klein 1AFC         ***')
z=erfinv(2*pHpF-1)*sqrt(2)      %converts p to z
DprimeKlein=-diff(z)            %dprime = zhit - zfalse alarm
CritC=-sum(z)/2                 %criterion= (zhit+zfalse alarm)/2
PcKlein= (pH + (1-pF))/2        %mean %correct(commonly used for 1AFC&2AFC)

%% Palamedes Method for 1AFC
disp('***      Palamedes 1AFC        ***')
[dp C lnBeta Pc]=PAL_SDT_1AFC_PHFtoDP(pHpF');
DprimePAL=dp'
CritC=C'
propCorrPal=Pc'

%% Palamedes Method for 2AFC
disp('***    Four 2AFC methods      ***')
[dP C lnB pC]=PAL_SDT_2AFC_PHFtoDP(pHpF'); %PAL method with bias allowed
DprimePAL=dP'
DprimeTimesSqrt2_ToCompareTo1AFC=dP'*sqrt(2) %to connect 1AFC
z2AFC=z/sqrt(2);            %Always that sqrt(2) to connect to 1AFC
DprimeKlein=-diff(z2AFC)    %same as Klein method above
pC_PAL=pC'                  %from the PAL 2AFC_PHFtoDP function
Dprime4pC = erfinv(2*pC'-1)*2        %simply double the z score above
dP_MAFC=PAL_SDT_MAFC_PCtoDP(pC',2)   %Using the MAFC routine for no bias
plot(x,[DprimeKlein;Dprime4pC],x,DprimeKlein,'.')
legend('including bias','ignoring bias','location','NorthWest')
title('dprime including and ignoring bias')
xlabel('stimulus strength');ylabel('dprime')
%Note that nowhere is Pc_max calculated. It is an artificial construct. 
disp('My convention is that there is only one pC, pCmax isn"t used')
disp('There are two types of dprime:  with and without bias')
