%This routine replicates any condition presented in Figure 3 of 'The
%Psychometric Function: The Lapse Rate Revisited' by Nick Prins.
%For each simulation it will display the likelihood function and the
%resulting best-fitting PF.
%At completion, it will create a figure akin to that shown in Figure3.pdf.
%This routine requires Palamedes version 1.4.1 or higher (available from
%www.palamedestoolbox.org)

clear all,

%s = RandStream.create('mt19937ar','seed',sum(100*clock)); %Do something different every time
%andStream.setDefaultStream(s);

warningstates = warning('query','all');%avoids inconsequential Log of Zero
warning off all                        %warnings                                         

doplot = 0;
s = 2
Ntrial = 120
numSims = 100
lapseGen = 0.01 
stepThru = 'n';

while Ntrial ~= 120 && Ntrial ~= 240 && Ntrial ~= 480 && Ntrial ~= 960
    Ntrial = input('N (120, 240, 480, or 960): ');
end
%numSims = input('Number of simulations (e.g., 2000): ');
while s < 1 || s > 7
    s = input('Placement regimen (1 through 7): ');
end
%lapseGen = input('Generating lapse rate (0 through 0.05): ');

%stepThru = input('Step through by spacebar (y/n): ','s');

switch Ntrial   %Vary histogram axis limits etc. with N
    case 120
        crit = log10([5 17; 3.5 22; .002 9; .0005 16]);        
    case 240
        crit = log10([6 14; 5 16; .008 2; .003 3]);        
    case 480
        crit = log10([7 12; 6 13; .02 1; .008 1.4]);        
    case 960
        crit = log10([7 11; 6.5 12; .03 .5; .02 .6]);        
end

%options for simplex search
options = PAL_minimize('options');
options.TolFun = 1e-09;     %low tolerance (i.e., high precision)
options.TolX = 1e-09;       %low tolerance (i.e., high precision)
options.MaxFunEvals = 10000;
options.MaxIter = 10000;
options.Display = 'off';

PF = @PAL_Weibull;          %Generating and fitted function

%Stimulus placement regimens. Values obtained from Hill, 2001 (Retrieved 
%from: http://bootstrap-software.org/psignifit/publications/hill2001.pdf):
s1 = PF([10 3 0 0],[.3 .4 .48 .52 .6 .7],'inverse');
s2 = PF([10 3 0 0],[.1 .3 .4 .6 .7 .9],'inverse');
s3 = PF([10 3 0 0],[.3 .44 .7 .8 .9 .98],'inverse');
s4 = PF([10 3 0 0],[.1 .2 .3 .4 .5 .6],'inverse');
s5 = PF([10 3 0 0],[.08 .18 .28 .7 .85 .99],'inverse');
s6 = PF([10 3 0 0],[.3 .4 .5 .6 .7 .99],'inverse');
s7 = PF([10 3 0 0],[.34 .44 .54 .8 .9 .98],'inverse');

OutOfNum = (Ntrial/6)*ones(1,6);        %Number of trials at each intensity

lapseLimits = [0 0.06];                 %prior on lambda

paramsFit = zeros(numSims,3,4);           %pre-allocate memory
paramsFitAdj = zeros(numSims,3,2);
LL = zeros(numSims,3);
exitflag = zeros(numSims,3);

%Placement

switch s
    case 1
        x = s1;
    case 2
        x = s2;
    case 3
        x = s3;
    case 4
        x = s4;
    case 5
        x = s5;
    case 6
        x = s6;
    case 7
        x = s7;
end

%for plotting later
xFine = PF([10 3 0 0],.05,'inverse'):(PF([10 3 0 0],.999,'inverse')-PF([10 3 0 0],.05,'inverse'))/100:PF([10 3 0 0],.999,'inverse');

%Create 3D parameter space (threshold, slope, lapse) and log probability
%look-up-tables to be used during brute-force search for seeds for simplex 
%search.
%Note: this could be performed within function PAL_PFML_Fit (or
%PAL_PFML_BruteForceFit) but is moved outside of simulation loop to reduce
%computational load. Type help PAL_PFML_Fit or help PAL_PFML_BruteForceFit 
%for more information.

gridGrain = 150;

alphas = logspace(log10(PAL_Weibull([10 3 0 0],.01,'inverse')),log10(PAL_Weibull([10 3 0 0],.99,'inverse')),gridGrain);
betas = logspace(-1,2,gridGrain);
gamma = 0.5;
lambdas = [0:.005:.06];

[params.alpha params.beta params.lambda] = ndgrid(alphas,betas,lambdas);
params.gamma = ones(size(params.alpha)).*gamma;

logpcorrect = zeros(length(x),size(params.alpha,1),size(params.alpha,2),size(params.alpha,3));
logpincorrect = zeros(length(x),size(params.alpha,1),size(params.alpha,2),size(params.alpha,3));
for level = 1:length(x)
    logpcorrect(level,:,:,:) = log(PF(params,x(1,level)));
    logpincorrect(level,:,:,:) = log(1-PF(params,x(1,level)));
end

paramsGen = [10 3 .5 lapseGen]; %generating parameters

if doplot==1 % drc
figure('units','pixels','position',[10 100 1300 500]);
end

[gridx gridy] = meshgrid(1:gridGrain,1:gridGrain);

%Simulation loop
for simulation = 1:numSims
    if mod(simulation,10) == 0    %for the impatient
        disp(int2str(simulation));
    end

    %simulate observer
    NumPos = PAL_PF_SimulateObserverParametric(paramsGen, x, OutOfNum, PF);

%%%%%%%%%%%%%%%%%Fit using free lapse rate

    %First perform brute force search for initial values to seed simplex
    %search with. This is moved outside of function PAL_PFML_Fit to reduce 
    %computational load (see above).
    
    %Create likelihood grid and find maximum in it.
    LLspace = zeros(size(params.alpha,1),size(params.alpha,2),size(params.alpha,3));
    for level = 1:6 %6 stimulus intensities
       LLspace = LLspace + NumPos(level).*squeeze(logpcorrect(level,:,:,:))+(OutOfNum(level)-NumPos(level)).*squeeze(logpincorrect(level,:,:,:));
    end    
    
    [trash I] = PAL_findMax(LLspace);
    paramsInit.alpha = alphas(I(1));
    paramsInit.beta = betas(I(2));
    paramsInit.gamma = .5;
    paramsInit.lambda = lambdas(I(3));    
    
    [LLspace2D Ilambda] = max(LLspace,[],3);
    
    LLspace2D = mat2gray(exp(LLspace2D));
        
    %fit function using simplex:
    paramsFree = [1 1 0 1];
    [paramsFit(simulation,1,:) LL(simulation,1) exitflag(simulation,1)] = PAL_PFML_Fit(x, NumPos, OutOfNum, paramsInit, paramsFree, PF,'lapseLimits',lapseLimits,'searchoptions',options);

    %%%%%%%Plot things
    
    if doplot==1 % drc
    %%%%%%%%%%%contour
    
    clf('reset')    
    contour(gridx,gridy,LLspace2D);
    set(gca,'units','pixels','position',[60 50 350 350]);
    xlabel('Slope (\beta)','fontsize',12);
    ylabel('Threshold (\alpha)','fontsize',12);
    zlabel('Likelihood','fontsize',12);

    axis([1 gridGrain 1 gridGrain]);

    line([I(2)-15 I(2)-15],[I(1)-15 I(1)+15],'color','k')
    line([I(2)+15 I(2)+15],[I(1)-15 I(1)+15],'color','k')
    line([I(2)-15 I(2)+15],[I(1)-15 I(1)-15],'color','k')
    line([I(2)-15 I(2)+15],[I(1)+15 I(1)+15],'color','k')

    S = [1 10 100];
    xtick = [1 ((log10(S)-log10(.1)).*((gridGrain-1)/(log10(100)-log10(.1))))];

    set(gca,'xtick',xtick,'xticklabel',{'.1','1','10','100'})

    T = [2.2 4 6 10 16];
    ytick = ((log10(T)-log10(PAL_Weibull([10 3 0 0],.01,'inverse'))).*((gridGrain-1)/(log10(PAL_Weibull([10 3 0 0],.99,'inverse'))-log10(PAL_Weibull([10 3 0 0],.01,'inverse')))));

    set(gca,'ytick',ytick,'yticklabel',{'2.2','4','6','10','16'})
    text(1,185,'Likelihood values across functions in Brute-Force Searchgrid:','fontsize',12)
    text(1,175,'Searchgrid contains 292,500 PFs: 150 \alphas x 150 \betas x 13 \lambdas','fontsize',12)
    text(1,165,'Likelihoods shown for PFs with that \lambda (of the 13) that gives highest likelihood','fontsize',8)
    text(1,155,'(Best fitting function in grid serves as initial guess for Simplex)','fontsize',8)
    set(gca,'fontsize',12);
    set(gca,'handlevisibility','off');

    %%%%%%%%%%%%%Close up

    LLspaceCloseUp = LLspace2D(max(1,I(1)-15):min(I(1)+15,gridGrain),max(1,I(2)-15):min(I(2)+15,gridGrain));

    surf(LLspaceCloseUp,Ilambda(max(1,I(1)-15):min(I(1)+15,gridGrain),max(1,I(2)-15):min(I(2)+15,gridGrain)))
    set(gca,'units','pixels','position',[500 20 400 400]);
    set(gca,'xtick',[],'ytick',[],'ztick',[]);
    set(gca,'zgrid','off','xgrid','off','ygrid','off');
    zlabel('Likelihood','fontsize',12);

    axis([1 31 1 31 0 1]);    
    text(-4,15,'Threshold (\alpha)','rotation',-29,'fontsize',12,'horizontalalignment','center');
    text(15,-4,'Slope (\beta)','rotation',14,'fontsize',12,'horizontalalignment','center');
    set(gca,'fontsize',12);

    set(gca,'handlevisibility','off');

    %%%%%%%%%%%%colorbar
    cb = repmat(1:14,[2 1]);
    pcolor(cb)
    set(gca,'xgrid','off','ygrid','off');
    set(gca,'xtick',[],'ytick',[]);
    for i = 0:12
        Tl = num2str(i/200,'%.3f');
        text(i+1.5,2.5,Tl,'rotation',45);
    end
    text(7,0,'Lapse rate','horizontalalignment','center','fontsize',12)
    text(7,5.5,'Close up of square in plot on left','horizontalalignment','center','fontsize',12)
    text(7,4.5,'(Color code indicates best lambda here)','horizontalalignment','center','fontsize',8)
    set(gca,'units','pixels','position',[500 400 400 20]);
    set(gca,'handlevisibility','off');

    %%%%%%%%PFs

    paramsFitFull = squeeze(paramsFit(simulation,1,:));

    semilogx(xFine,PF(paramsGen,xFine),'-','color',[0 0 0],'linewidth',2)
    hold on
    semilogx(xFine,PF(paramsInit,xFine),'-','color',[0 0 1],'linewidth',3)
    semilogx(xFine,PF(paramsFitFull,xFine),'-','color',[1 0 0],'linewidth',2)

    axis([3 20 .4 1])   
    semilogx(x,NumPos/(Ntrial/6),'o','color',[0 0 0],'markersize',8,'markerfacecolor','k')
    ylabel('Proportion Correct','fontsize',16);
    xlabel('Stimulus Intensity','fontsize',16);
    set(gca,'xtick',[4:19],'xticklabel',{'4','5','6','7','8','9','','11','','13','','','16','','','19'});
    set(gca,'ytick',[.5:.1:1]);
    set(gca,'units','pixels','position',[975 75 300 300])

    text(min(xFine),1.2,'Fits:','fontsize',12)
    text(min(xFine),1.15,'Generating function','fontsize',12)
    text(min(xFine),1.1,'Best function in discrete searchgrid','fontsize',12,'color',[0 0 1])
    text(min(xFine),1.05,'Best function by Simplex','fontsize',12,'color',[1 0 0])
    
    set(gca,'handlevisibility','off')

    semilogx(xFine,PF([10 3 0 0],xFine),'-','color',[0 0 0],'linewidth',1)
    axis([3 20 0 1])
    hold on
    semilogx(xFine,PF([paramsInit.alpha paramsInit.beta 0 0],xFine),'-','color',[0 0 1],'linewidth',2)
    semilogx(xFine,PF([paramsFitFull(1) paramsFitFull(2) 0 0],xFine),'-','color',[1 0 0],'linewidth',1)
    set(gca,'units','pixels','position',[1000 250 100 100])
    set(gca,'ytick',[0:1]);
    set(gca,'xtick',[4:19],'xticklabel',{'4','','','','','','','','','','','','','','','19'});
    text(2.5,.5,'F(\itx\rm)','horizontalalignment','center','rotation',90)
    set(gca,'handlevisibility','off')

    drawnow

    if strcmpi(stepThru,'y')
        pause
    end

end % drc / doplot
    
    %transform to metric reported by Wichmann & Hill
    paramsFitAdj(simulation,1,1) = PF([paramsFit(simulation,1,1) paramsFit(simulation,1,2) 0 0],.5,'inverse');
    paramsFitAdj(simulation,1,2) = PF([paramsFit(simulation,1,1) paramsFit(simulation,1,2) 0 0],paramsFitAdj(simulation,1,1),'derivative');

%%%%%%%%%%%%%%%Fit using fixed lapse rate = 0                                        

    lapseFit = 0;
    LLspace = zeros(size(params.alpha,1),size(params.alpha,2));
    for level = 1:6
       LLspace = LLspace + NumPos(level).*squeeze(logpcorrect(level,:,:,find(lambdas==lapseFit)))+(OutOfNum(level)-NumPos(level)).*squeeze(logpincorrect(level,:,:,find(lambdas==lapseFit)));
    end
    
    [trash I] = PAL_findMax(LLspace);
    paramsInit.alpha = alphas(I(1));
    paramsInit.beta = betas(I(2));
    paramsInit.gamma = .5;
    paramsInit.lambda = lapseFit;
    
    %fit function:
    paramsFree = [1 1 0 0];
    [paramsFit(simulation,2,:) LL(simulation,2) exitflag(simulation,2)] = PAL_PFML_Fit(x, NumPos, OutOfNum, paramsInit, paramsFree, PF,'searchoptions',options);
    
    %transform to metric reported by Wichmann & Hill
    paramsFitAdj(simulation,2,1) = PF([paramsFit(simulation,2,1) paramsFit(simulation,2,2) 0 0],.5,'inverse');
    paramsFitAdj(simulation,2,2) = PF([paramsFit(simulation,2,1) paramsFit(simulation,2,2) 0 0],paramsFitAdj(simulation,2,1),'derivative');

%%%%%%%%%%%%%%%Fit using fixed lapse rate = 0.025

    lapseFit = 0.025;
    LLspace = zeros(size(params.alpha,1),size(params.alpha,2));
    for level = 1:6
       LLspace = LLspace + NumPos(level).*squeeze(logpcorrect(level,:,:,find(lambdas==lapseFit)))+(OutOfNum(level)-NumPos(level)).*squeeze(logpincorrect(level,:,:,find(lambdas==lapseFit)));
    end
    [trash I] = PAL_findMax(LLspace);
    paramsInit.alpha = alphas(I(1));
    paramsInit.beta = betas(I(2));
    paramsInit.gamma = .5;
    paramsInit.lambda = lapseFit;

    %fit function:
    paramsFree = [1 1 0 0];
    [paramsFit(simulation,3,:) LL(simulation,3) exitflag(simulation,3)] = PAL_PFML_Fit(x, NumPos, OutOfNum, paramsInit, paramsFree, PF,'searchoptions',options);

    %transform to metric reported by Wichmann & Hill
    paramsFitAdj(simulation,3,1) = PF([paramsFit(simulation,3,1) paramsFit(simulation,3,2) 0 0],.5,'inverse');
    paramsFitAdj(simulation,3,2) = PF([paramsFit(simulation,3,1) paramsFit(simulation,3,2) 0 0],paramsFitAdj(simulation,3,1),'derivative');

end

%Display results
%% return warningstates to original settings
warning(warningstates);
name = ['Figure 3, N = ', int2str(Ntrial), ', Placement: s', int2str(s), ', Generating lapse rate: ', num2str(lapseGen)];
figure('units','pixels','position',[100 10 750 565],'name',name);

nonpos = paramsFitAdj <= 0;
paramsFitAdj = log10(paramsFitAdj);
paramsFitAdj(nonpos) = -Inf;

for cond = 1:3
    paramsFitAdj(:,cond,3:4) = paramsFit(:,cond,3:4);
    paramsFitAdj(paramsFitAdj(:,cond,1)>=crit(1,2),cond,1) = crit(2,2);
    paramsFitAdj(paramsFitAdj(:,cond,1)<= crit(1,1),cond,1) = crit(2,1);
    paramsFitAdj(paramsFitAdj(:,cond,2)>=crit(3,2),cond,2) = crit(4,2);
    paramsFitAdj(paramsFitAdj(:,cond,2)<= crit(3,1),cond,2) = crit(4,1);
end

lolapse = paramsFitAdj(paramsFitAdj(:,1,4)<.0001,1,:);
hilapse = paramsFitAdj(paramsFitAdj(:,1,4)>.0599,1,:);
melapse = paramsFitAdj(paramsFitAdj(:,1,4)>.0001&paramsFitAdj(:,1,4)<.0599,1,:);

%Lapse Rates

lapseBinCenters = .06/80:.06/40:.06-.06/80;
hist(squeeze(paramsFitAdj(:,1,4)),lapseBinCenters);
set(gca,'Units','pixels','position',[53 473 289 49],'fontsize',8);
maxim = max(hist(squeeze(paramsFitAdj(:,1,4)),lapseBinCenters));
XL = [0 0.06];
axis([XL 0 maxim*1.3]);
set(gca,'xtick',[0:.01:.06]);
message = ['N = ' int2str(size(paramsFitAdj(:,1,4),1))];
text(.0025,maxim,message,'fontsize',8);
xlabel('lapse rate estimate','fontsize',8);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[.5 .5 .5])
line(XL,[0 0],'linewidth',1,'color','k');
set(gca,'handlevisibility','off');
plot(lapseGen,maxim*1.22,'kv','markerfacecolor','k');
set(gca,'Units','pixels','position',[53 473 289 49]);
axis([XL 0 maxim*1.3]);
text(mean(XL),maxim*1.5,'\bf lapse rate free to vary','horizontalalignment','center','fontsize',8);        
axis off
set(gca,'handlevisibility','off');

%%%%%%%%%%%%%%%Thresholds

binCenters = log10(logspace(crit(2,1),(crit(2,2)),40));

for histogram = 1:6

    switch histogram
        case 1
            toPlot = squeeze(paramsFitAdj(:,1,1));
            toPlotMessage = ['\bfAll'];
            condMessage = [''];
            yPos = 390;
        case 2
            toPlot = squeeze(lolapse(:,1));
            toPlotMessage = ['\bf\lambda_{est} = 0'];
            condMessage = [''];
            yPos = 338;
        case 3
            toPlot = squeeze(hilapse(:,1));
            toPlotMessage = ['\bf\lambda_{est} = 0.06'];
            condMessage = [''];
            yPos = 285;
        case 4
            toPlot = squeeze(melapse(:,1));
            toPlotMessage = ['\bfother'];
            condMessage = [''];
            yPos = 233;
        case 5
            toPlot = squeeze(paramsFitAdj(:,2,1));
            toPlotMessage = ['\bfAll'];
            condMessage = ['\bflapse rate fixed at 0'];
            yPos = 139;
        case 6
            toPlot = squeeze(paramsFitAdj(:,3,1));
            toPlotMessage = ['\bfAll'];
            condMessage = ['\bflapse rate fixed at 0.025'];
            yPos = 45;
    end

    if ~isempty(toPlot)
        hist(toPlot,binCenters);
        set(gca,'Units','pixels','position',[53 yPos 289 49],'fontsize',8);
        maxim = max(hist(toPlot,binCenters));
        XL = [crit(2,1)-(crit(2,2)-crit(2,1))/20 crit(2,2)+(crit(2,2)-crit(2,1))/20];
        axis([XL 0 maxim*1.3]);
        xtick = get(gca,'xtick');
        xtick(1) = crit(2,1);
        xtick(length(xtick)) = crit(2,2);
        xtick(2:length(xtick)-1) = log10(round(10.^xtick(2:length(xtick)-1)));
        set(gca,'xtick',xtick);
        xticklabel = {};
        for index = 1:length(xtick); xticklabel(index) = cellstr(num2str(10.^xtick(index)));end
        xticklabel(1) = cellstr(['<' num2str(10.^crit(1,1))]);
        xticklabel(length(xticklabel)) = cellstr(['>' num2str(10.^crit(1,2))]);
        set(gca,'xticklabel',xticklabel);
        text(xtick(1),maxim*1.1,toPlotMessage,'fontsize',8);
        message = ['N = ' int2str(size(toPlot,1))];
        text(xtick(1),maxim*.75,message,'fontsize',8);
        Median = 10.^median(toPlot);
        message = ['Median = ' num2str(Median,'%5.3f')];
        text(xtick(length(xtick)),maxim,message,'horizontalalignment','right','fontsize',8);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[.5 .5 .5])
        line(XL,[0 0],'linewidth',1,'color','k');
        xlabel('Threshold estimate','fontsize',8);
        text(mean(XL),maxim*1.5,condMessage,'horizontalalignment','center','fontsize',8);
        set(gca,'handlevisibility','off');
        plot(log10(8.85),maxim*1.22,'kv','markerfacecolor','k');
        set(gca,'Units','pixels','position',[53 yPos 289 49]);
        axis([XL 0 maxim*1.3]);
        axis off
        set(gca,'handlevisibility','off');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%Slopes

binCenters = log10(logspace(crit(4,1),crit(4,2),40));

for histogram = 1:6

    switch histogram
        case 1
            toPlot = squeeze(paramsFitAdj(:,1,2));
            toPlotMessage = ['\bfAll'];
            condMessage = [''];
            yPos = 390;
        case 2
            toPlot = squeeze(lolapse(:,2));
            toPlotMessage = ['\bf\lambda_{est} = 0'];
            condMessage = [''];
            yPos = 338;
        case 3
            toPlot = squeeze(hilapse(:,2));
            toPlotMessage = ['\bf\lambda_{est} = 0.06'];
            condMessage = [''];
            yPos = 285;
        case 4
            toPlot = squeeze(melapse(:,2));
            toPlotMessage = ['\bfother'];
            condMessage = [''];
            yPos = 233;
        case 5
            toPlot = squeeze(paramsFitAdj(:,2,2));
            toPlotMessage = ['\bfAll'];
            condMessage = ['\bflapse rate fixed at 0'];
            yPos = 139;
        case 6
            toPlot = squeeze(paramsFitAdj(:,3,2));
            toPlotMessage = ['\bfAll'];
            condMessage = ['\bflapse rate fixed at 0.025'];
            yPos = 45;
    end
    
    if ~isempty(toPlot)
    
        hist(toPlot,binCenters);
        set(gca,'Units','pixels','position',[428 yPos 289 49],'fontsize',8);
        maxim = max(hist(toPlot,binCenters));
        XL = [crit(4,1)-(crit(4,2)-crit(4,1))/20 crit(4,2)+(crit(4,2)-crit(4,1))/20];
        axis([XL 0 maxim*1.3]);
        xtick = get(gca,'xtick');
        xtick = xtick(10.^xtick >= .01);
        xtick(1) = crit(4,1);
        xtick(length(xtick)) = crit(4,2);
        xtick(2:length(xtick)-1) = log10((round(100*(10.^xtick(2:length(xtick)-1))))/100);
        set(gca,'xtick',xtick);
        xticklabel = {};
        for index = 1:length(xtick); xticklabel(index) = cellstr(num2str(10.^xtick(index)));end
        xticklabel(1) = cellstr(['<' num2str(10.^crit(3,1))]);
        xticklabel(length(xticklabel)) = cellstr(['>' num2str(10.^crit(3,2))]);
        set(gca,'xticklabel',xticklabel);
        text(xtick(1),maxim*1.1,toPlotMessage,'fontsize',8);
        message = ['N = ' int2str(size(toPlot,1))];
        text(xtick(1),maxim*.75,message,'fontsize',8);
        Median = 10.^median(toPlot);
        message = ['Median = ' num2str(Median,'%5.3f')];
        text(xtick(length(xtick)),maxim,message,'horizontalalignment','right','fontsize',8);
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[.5 .5 .5])
        line(XL,[0 0],'linewidth',1,'color','k');
        xlabel('Slope estimate','fontsize',8);
        text(mean(XL),maxim*1.5,condMessage,'horizontalalignment','center','fontsize',8);
        set(gca,'handlevisibility','off');
        plot(log10(.118),maxim*1.22,'kv','markerfacecolor','k');
        set(gca,'Units','pixels','position',[428 yPos 289 49]);
        axis([XL 0 maxim*1.3]);
        axis off
        set(gca,'handlevisibility','off');
    end
    
end

%% get SE and correlations
params=squeeze(paramsFit(:,1,[1 2 4]));
meanParams=mean(params)
covar=cov(params);
SE=sqrt(diag(covar))'
cor=covar./(SE'*SE);
correl=[cor(1,2) cor(1,3) cor(2,3)]

paramsAdj=10.^squeeze(paramsFitAdj(:,1,[1 2 4]));
meanParamsAdj=mean(paramsAdj)
covar=cov(paramsAdj);
SEAdj=sqrt(diag(covar))'
cor=covar./(SEAdj'*SEAdj);
correlAdj=[cor(1,2) cor(1,3) cor(2,3)]
