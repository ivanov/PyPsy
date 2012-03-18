clear all;  format compact
Homework=3,   %1) Ex 2a, p. 222, 2) p 225, Ex2a,  4) 2b, 5) 3a, 6) 3b
B=200;   %this will give p values in steps of 0.005 
tic
PF = @PAL_Logistic;
if Homework<5
    StimLevels = [-2:1:2; -2:1:2];
    OutOfNum = [100 100 100 100 100; 100 100 100 100 100];
    NumPos = [61 70 81 92 97; 59 59 67 86 91];
    PF = @PAL_Logistic;
    params = [0 1 .5 0; 0 1 .5 0];
    if Homework==1,  %same as p. 222 
        [TLR pTLR paramsL paramsF TLRSim converged] = PAL_PFLR_ModelComparison ...
            (StimLevels, NumPos, OutOfNum, params, B, PF);
        paramsL
        paramsF
        pTLR    %probability of transformed likelihood ratio
        
    elseif Homework==2  % same as p. 225
        %exercise 2a again  but this time uncontrain the lesser
        [TLR pTLR paramsL paramsF TLRSim converged] = PAL_PFLR_ModelComparison ...
            (StimLevels, NumPos, OutOfNum, params, B, PF, 'lesserSlopes', 'unconstrained');
        paramsL
        paramsF
        pTLR
    elseif Homework==3  %exercise 2a but now constrain fuller slopes
        [TLR pTLR paramsL paramsF TLRSim converged] = PAL_PFLR_ModelComparison ...
            (StimLevels, NumPos, OutOfNum, params, B, PF, 'fullerSlopes', 'constrained');
        paramsL
        paramsF
        pTLR
    elseif Homework==4
        %exercise 2b
        [TLR pTLR paramsL paramsF TLRSim converged] = PAL_PFLR_ModelComparison ...
            (StimLevels, NumPos, OutOfNum, params, B, PF, 'fullerThresholds', 'constrained');
        paramsL
        paramsF
        pTLR
    end
end

if Homework==5,
    %Exercise 3a    contrast: 0 seconds V.S. 4 seconds
    StimLevels = [-2:1:2; -2:1:2];
    OutOfNum = [150 150 150 150 150; 150 150 150 150 150];
    NumPos = [84 92 128 137 143; 67 85 103 131 139];
    PF = @PAL_Logistic;
    params = [-.6 1.8 .5 .02; .1 1.8 .5 .02];
    [TLR pTLR paramsL paramsF TLRSim converged] = ...
        PAL_PFLR_ModelComparison(StimLevels, NumPos, OutOfNum, ...
        params, B, PF, 'lesserLapseRates', 'constrained', 'fullerLapseRates', ...
        'constrained', 'fullerSlopes', 'constrained', 'lapseLimits',[0 1]);
    pTLR
end

if Homework==6,  %contrast: 4 second V.S. 8 seconds V.S. 12 seconds
    StimLevels = [-2:1:2; -2:1:2; -2:1:2];
    OutOfNum = [150 150 150 150 150; 150 150 150 150 150; 150 150 150 150 150];
    NumPos = [67 85 103 131 139; 73 85 92 125 143; 82 86 97 122 141];
    params = [.1 1.8 .5 .02; .6 1.8 .5 .02; .9 1.8 .5 .02];
    
    [paramsFitted LL exitflag output] = ...
        PAL_PFML_FitMultiple(StimLevels, NumPos, OutOfNum, params, PF, 'Thresholds', 'unconstrained', ...)
        'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates', 'constrained');
    paramsFitted
    [TLR pTLR TLRSim coveraged] = ...
        PAL_PFML_GoodnessOfFitMultiple(StimLevels, NumPos, OutOfNum, paramsFitted, B, PF, 'Thresholds', ...
        'unconstrained', 'Slopes', 'constrained', 'GuessRates', 'fixed', 'LapseRates', 'constrained', 'lapseLimits', [0 1]);
    pTLR
end
toc

