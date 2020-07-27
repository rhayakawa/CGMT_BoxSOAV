%%--------------------------------------------------------------------------------------------
% comparison with matched filter bound
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
Delta=1;
p0=0.9;
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
L=length(arrR);
arrSNR=0:1:20;
q2Iteration=50;

arrSER_theo_naive=zeros(size(arrSNR));
arrSER_theo_optimal=zeros(size(arrSNR));
arrSER_MFB=zeros(size(arrSNR));
for SNRIndex=1:length(arrSNR)
  SNR=arrSNR(SNRIndex);
  disp(['SNR=' num2str(SNR)]);
  sigma2_v=arrP*(arrR.^(2)).'/(10^(SNR/10));  % noise variance

  % get optimal parameters and quantizer
  q2Min=0;
  q2Max=1;
  for q2IterationIndex=1:q2Iteration
    disp(['  q2=' num2str(q2Min) '-' num2str(q2Max)]);
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
  
  % naive quantizer with optimal q2
  q2Min=0;
  q2Max=1;
  SER1=0;
  SER2=Inf;
  arrThr=[-Inf (arrR(1)+arrR(2))/2 (arrR(2)+arrR(3))/2 Inf];
  for q2IterationIndex=1:q2Iteration
    disp(['  q2=' num2str(q2Min) '-' num2str(q2Max)]);
    q2_1=(2*q2Min+q2Max)/3;
    q2_2=(q2Min+2*q2Max)/3;
    % SER for q2_1
    arrCoef=[1 q2_1 1];
    arrQ=arrCoef*[-1  1  1  1;...
                  -1 -1  1  1;...
                  -1 -1 -1  1;];
    arrQ(1)=-Inf;
    arrQ(L+1)=Inf;
    [~,SER1,~,~]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
    % SER for q2_2
    arrCoef=[1 q2_2 1];
    arrQ=arrCoef*[-1  1  1  1;...
                  -1 -1  1  1;...
                  -1 -1 -1  1;];
    arrQ(1)=-Inf;
    arrQ(L+1)=Inf;
    [~,SER2,~,~]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
    if SER1<SER2
      q2Max=q2_2;
    else
      q2Min=q2_1;
    end
    if abs(q2_2-q2_1)/abs(q2_2)<1e-4
      break;
    end
  end
  arrSER_theo_naive(SNRIndex)=SER1;

  % matched filter bound
  kappa2_MFB=sigma2_v*log((1-p0)/(2*p0))-Delta/2;
  arrSER_MFB(SNRIndex)...
    =1-((1-p0)*normcdf(1/(sqrt(Delta*sigma2_v))*(Delta+kappa2_MFB))...
        +p0*(normcdf(-kappa2_MFB/sqrt(Delta*sigma2_v))-normcdf(kappa2_MFB/sqrt(Delta*sigma2_v))));
end

%% plot result
close all;
figure;
setLegend={};
h=semilogy(arrSNR,arrSER_theo_naive,'--k','LineWidth',1);
hold on;
setLegend=[setLegend '$\mathcal{Q}_{\mathrm{NV}}(\cdot)$'];
h=semilogy(arrSNR,arrSER_theo_optimal,'-k','LineWidth',1);
setLegend=[setLegend '$\mathcal{Q}_{\mathrm{AO}}(\cdot)$'];
h=semilogy(arrSNR,arrSER_MFB,':k','LineWidth',1);
setLegend=[setLegend 'MFB'];
grid on;

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
ylabel('asymptotic SER');
fig.XLim=[min(arrSNR) max(arrSNR)];
fig.YLim=[1e-6 1e-1];
fig.XTick=0:5:30;
saveas(h, ['SERvsSNR_TS_MFB(p0=' num2str(p0) ',Delta=' num2str(Delta) ').eps'], 'epsc');
