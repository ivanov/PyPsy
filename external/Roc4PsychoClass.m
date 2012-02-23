function ROC4PsychoClass(data, params0, info)
disp('****************Rating STD analysis*******************')
LLorChisq=info(3);  %2: for chi square minimumization 
[Nresp Nstim]=size(data)
CumSum=cumsum(data);
Total=CumSum(end,:); %total number of trials at each level
CumProb=1-CumSum./(ones(Nresp,1)*Total);%this make probs to RIGHT of criterion
CumProb=[ones(1,Nstim); CumProb] %the extra row is for pretty plot
subplot(1,2,1)
plot(CumProb(:,1),CumProb(:,2:4),'-',CumProb(:,1),CumProb(:,2:4),'*', ...
    [0 1],[0 1], '-');
legend('low','medium','high','location','southeast');hold on
z=p2zfun(CumProb)
%z=erfinv(2*CumProb-1)*sqrt(2)  %each z row is for one criterion5
dprime=z-z(:,1)*ones(1,Nstim)  %recall the 1st column is for blanks
subplot(1,2,2)
plot(z(:,1), z(:,2:4),'-',z(:,1), z(:,2:4),'*',[-2 0],[-2 0])
xlabel('z score for false alarm (blank stimuli)')
ylabel('z score for hit rate for signals')
title('ROC curve for in z-score units')
legend('low','medium','high','location','southeast')

if Nresp==2 
    prob=data(2,:)./CumSum(end,:)
    z=erfinv(2*prob-1)*sqrt(2)
    dprime2resp=diff(z)
    disp('***   Palamedes 1AFC but note order of HF is reversed ***')
    [dp C lnBeta Pc]=PAL_SDT_1AFC_PHFtoDP(prob(1:2));
    DprimePAL=dp'
end

if LLorChisq==1,
    disp('this is LogLikelihood')
    paramsSimple = fminsearch('SetRocSimple',params0,[], offset, data,info)
    [LL EstSimpleLL]=SetRocSimple(paramsSimple, offset, data,info);
end
[paramsSimple,chisq1,f,EXITFLAG,OUTPUT,LAMBDA,j] = lsqnonlin('SetRocSimple',params0,[],[],[], data,info);
disp('Minimizing chisq')
chisq1
paramsSimple
j=full(j);
cov=inv(j'*j);
SE=sqrt(diag(cov))'
correl=cov./(SE'*SE)
[err, expectSimple]=SetRocSimple(paramsSimple, data,info);
expect=num2str(expectSimple,3)  
chisq=sum(sum((expectSimple-data).^2./expectSimple))
DegFree=(Nresp-1)*Nstim-length(params0)
subplot(1,2,1);
ProbPlotx=[.01:.01:1]; 
zPlot=p2zfun(ProbPlotx);
for i=1:Nstim-1;
    ProbPloty=z2pfun(zPlot+paramsSimple(i));
    plot(ProbPlotx,ProbPloty); hold on
end
xlabel('false alarm probability (blank stimuli)')
ylabel('hit rate for signals')
title('ROC curve for probability')


function z=p2zfun(p)
z=erfinv(2*p-1)*sqrt(2);

function  p = z2pfun(z)
p = [erf(z/sqrt(2)) + 1]/2;
    



