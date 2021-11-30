%%--------------------------------------------------------------------------------------------
% asymptotic SER versus SNR in discrete-valued sparse vector reconstruction
%
% Author:
%   Ryo Hayakawa
% Article:
%   Ryo Hayakawa and Kazunori Hayashi,
%   "Asymptotic performance of discrete-valued vector reconstruction 
%    via box-constrained optimization with sum of l1 regularizers,"
%   IEEE Transactions on Signal Processing, vol. 68, pp. 4320-4335, 2020. 
%%--------------------------------------------------------------------------------------------


clear;
addpath('subfunctions');

% parameters
N=1000;
Delta=0.8;
p0=0.8;
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
L=length(arrR);
arrSNR=0:2.5:30;
nSample=100;
nIteration=300;
q2Iteration=50;

arrq2_original=[0.01 0.1];

rng('shuffle');

arrSERvsSNR_original=zeros(length(arrq2_original),length(arrSNR));
arrSERvsSNR_optimal=zeros(size(arrSNR));
arrSER_theo_original=zeros(length(arrq2_original),length(arrSNR));
arrSER_theo_optimal=zeros(size(arrSNR));
for SNRIndex=1:length(arrSNR)
  SNR=arrSNR(SNRIndex);
  disp(['SNR=' num2str(SNR)]);

  % get optimal parameters and quantizer
  q2Min=0;
  q2Max=1;
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
    if abs(q2_2-q2_1)/abs(q2_2)<1e-4
      break;
    end
  end
  arrSER_theo_optimal(SNRIndex)=SER1;
  q2_opt=q2_1;
  arrCoef_opt=[1 q2_opt 1];
  arrQ_opt=arrCoef_opt*[-1  1  1  1;...
                        -1 -1  1  1;...
                        -1 -1 -1  1;];
  arrQ_opt(1)=-Inf;
  arrQ_opt(L+1)=Inf;
  [~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
  kappa2_opt=softThr(-1/2+alpha_opt^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
  arrThr_opt=[-Inf kappa2_opt -kappa2_opt Inf];
  
  [~,arrSERvsSNR_optimal(SNRIndex)]=get_empirical(N,Delta,arrP,arrR,arrQ_opt,arrThr_opt,SNR,nIteration,nSample);
  
  % original (Box-)SOAV
  for q2index=1:length(arrq2_original)
    q2=arrq2_original(q2index);
    arrCoef_original=[1 q2 1];
    arrQ_original=arrCoef_original*[-1  1  1  1;...
                                    -1 -1  1  1;...
                                    -1 -1 -1  1;];
    arrQ_original(1)=-Inf;
    arrQ_original(L+1)=Inf;
    arrThr_original=[-Inf -0.5 0.5 Inf];
    [~,arrSER_theo_original(q2index,SNRIndex),~,~]=get_theoretical(Delta,arrP,arrR,arrCoef_original,arrQ_original,SNR);
    [~,arrSERvsSNR_original(q2index,SNRIndex)]=get_empirical(N,Delta,arrP,arrR,arrQ_original,arrThr_original,SNR,nIteration,nSample);
  end
end

%% plot result
close all;
figure;
arrMarkers=['xk';'+k';'^k'];
arrLines=['--k';':k ';'-.k'];
setLegend={};
for q2index=1:length(arrq2_original)
  h=semilogy(arrSNR,arrSERvsSNR_original(q2index,:),arrMarkers(q2index,:),'LineWidth',1,'MarkerSize',10);
  setLegend=[setLegend ['$q_{2}=' num2str(arrq2_original(q2index)) '$ (empirical)']];
  hold on;
  h=semilogy(arrSNR,arrSER_theo_original(q2index,:),arrLines(q2index,:),'LineWidth',1);
  setLegend=[setLegend ['$q_{2}=' num2str(arrq2_original(q2index)) '$ (theoretical)']];
end
h=semilogy(arrSNR,arrSERvsSNR_optimal,'ok','LineWidth',1,'MarkerSize',10);
h=semilogy(arrSNR,arrSER_theo_optimal,'-k','LineWidth',1);
grid on;

setLegend=[setLegend 'optimal (empirical)' 'optimal (theoretical)'];
objLegend=legend(setLegend,'Location','northeast');
objLegend.Interpreter='latex';
objLegend.FontSize=16;
legend('boxoff');

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('SNR (dB)');
ylabel('SER');
fig.XLim=[min(arrSNR) max(arrSNR)];
fig.YLim=[1e-5 1e0];
fig.XTick=0:5:30;
saveas(h, ['SERvsSNR_TS(N=' num2str(N) ',p0=' num2str(p0) ',Delta=' num2str(Delta)  ',nSample=' num2str(nSample) ').eps'], 'epsc');
