%%--------------------------------------------------------------------------------------------
% empirical SER versus measurement ratio in discrete-valued vector reconstruction
%
% Author:
%   Ryo Hayakawa
% Article:
%   Ryo Hayakawa and Kazunori Hayashi,
%   "Asymptotic performance of discrete-valued vector reconstruction 
%    via box-constrained optimization with sum of l1 regularizers,"
%   IEEE Transactions on Signal Processing, vol. XX, no. XX, pp. XX-XX, 2020. 
%%--------------------------------------------------------------------------------------------


clear;
addpath('subfunctions');

% parameters
N=1000;
arrDelta=0.7:0.05:01;
p0=0.5;
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
L=length(arrR);
arrThr_naive=[-Inf -0.5 0.5 Inf];
SNR=20;
q2Iteration=50;
nIteration=200;
nSample=100;

arrSER_empi=zeros(1,length(arrDelta));
arrSER_SOAV_empi=zeros(1,length(arrDelta));
arrSER_BOX_empi=zeros(1,length(arrDelta));
arrSER_theo=zeros(1,length(arrDelta));
arrq2opt=zeros(1,length(arrDelta));
for DeltaIndex=1:length(arrDelta)
  Delta=arrDelta(DeltaIndex);
  disp(['Delta=' num2str(Delta)]);

  q2Min=0;
  q2Max=0.1;
  SER1=0;
  SER2=Inf;
  for q2IterationIndex=1:q2Iteration
    q2_1=(2*q2Min+q2Max)/3;
    q2_2=(q2Min+2*q2Max)/3;
    % SER for q2_1
    arrCoef=[1 q2_1 1];
    arrQ=arrCoef*[-1  1  1  1;...
                  -1 -1  1  1;...
                  -1 -1 -1  1;];
    arrQ(1)=-Inf;
    arrQ(L+1)=Inf;
    [~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
    kappa2=softThr(-1/2+alpha_opt^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
    arrThr=[-Inf kappa2 -kappa2 Inf];
    SER1=SER_theo(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ,arrThr);
    % SER for q2_2
    arrCoef=[1 q2_2 1];
    arrQ=arrCoef*[-1  1  1  1;...
                  -1 -1  1  1;...
                  -1 -1 -1  1;];
    arrQ(1)=-Inf;
    arrQ(L+1)=Inf;
    [~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
    kappa2=softThr(-1/2+alpha_opt^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
    arrThr=[-Inf kappa2 -kappa2 Inf];
    SER2=SER_theo(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ,arrThr);
    if SER1<SER2
      q2Max=q2_2;
    else
      q2Min=q2_1;
    end
    if abs(SER2-SER1)/abs(SER2)<1e-4
      break;
    end
  end
  arrSER_theo(DeltaIndex)=SER1;

  q2_opt=q2Min;
  arrq2opt(DeltaIndex)=q2_opt;
  arrCoef=[1 q2_opt 1];
  arrQ=arrCoef*[-1  1  1  1;...
                -1 -1  1  1;...
                -1 -1 -1  1;];
  [~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
  % threshold of AOQ
  kappa2=softThr(-1/2+alpha_opt^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
  arrThr=[-Inf kappa2 -kappa2 Inf];

  [MSE,SER]=get_empirical(N,Delta,arrP,arrR,arrQ,arrThr,SNR,nIteration,nSample);
  [MSE_SOAV,SER_SOAV]=get_empirical_noBOX(N,Delta,arrP,arrR,arrQ,arrThr,SNR,nIteration,nSample);
  [MSE_BOX,SER_BOX]=get_empirical(N,Delta,arrP,arrR,zeros(1,4),arrThr_naive,SNR,nIteration,nSample);

  arrSER_empi(DeltaIndex)=SER;
  arrSER_SOAV_empi(DeltaIndex)=SER_SOAV;
  arrSER_BOX_empi(DeltaIndex)=SER_BOX;
  
end

%% plot SER
close all;
figure;
arrLegend={};
% empirical
h=semilogy(arrDelta,arrSER_BOX_empi,'sk','LineWidth',1,'MarkerSize',10);
arrLegend=[arrLegend ['Box relaxation (empirical)']];
hold on;
h=semilogy(arrDelta,arrSER_SOAV_empi,'xk','LineWidth',1,'MarkerSize',10);
arrLegend=[arrLegend ['SOAV (empirical)']];
h=semilogy(arrDelta,arrSER_empi,'^k','LineWidth',1,'MarkerSize',10);
arrLegend=[arrLegend ['Box-SOAV (empirical)']];
grid on;

xlabel('$\Delta$');
ylabel('SER');
objLegend=legend(arrLegend);
objLegend.Interpreter='latex';
objLegend.Location='southwest';
objLegend.FontSize=20;
legend('boxoff');
fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
axis([min(arrDelta) max(arrDelta) 1e-5 1]);
saveas(h, ['SERvsDelta_TS(N=' num2str(N) ',p0=' num2str(p0) ',SNR=' num2str(SNR) ').eps'], 'epsc');
