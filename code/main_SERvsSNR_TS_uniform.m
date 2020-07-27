%%--------------------------------------------------------------------------------------------
% asymptotic SER for uniformly distributed vector
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
Delta=0.7;
p0=1/3;
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
L=length(arrR);
arrSNR=0:30;
q2Iteration=50;

arrSER_theo_uniform=zeros(size(arrSNR));
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
  [~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ_opt,SNR);
  kappa2_opt=softThr(-1/2+alpha_opt^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ_opt,arrR);
  arrThr_opt=[-Inf kappa2_opt -kappa2_opt Inf];
  
  % original (Box-)SOAV
  arrCoef_original=[1 0 1];
  arrQ_original=arrCoef_original*[-1  1  1  1;...
                                  -1 -1  1  1;...
                                  -1 -1 -1  1;];
  arrQ_original(1)=-Inf;
  arrQ_original(L+1)=Inf;
  arrThr_original=[-Inf -0.5 0.5 Inf];
  [~,arrSER_theo_uniform(SNRIndex),~,~]=get_theoretical(Delta,arrP,arrR,arrCoef_original,arrQ_original,SNR);
end

%% plot result
close all;
figure;
h=semilogy(arrSNR,arrSER_theo_uniform,'--xk','LineWidth',1,'MarkerSize',12,'MarkerIndices',1:5:31);
hold on;
h=semilogy(arrSNR,arrSER_theo_optimal,'-k','LineWidth',1);
grid on;

setLegend={'$q_{2}=0$','optimal $q_{2}$'};
objLegend=legend(setLegend,'Location','northeast');
objLegend.Interpreter='latex';
objLegend.FontSize=18;
legend('boxoff');

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('SNR (dB)');
ylabel('asymptotic SER');
fig.XLim=[min(arrSNR) max(arrSNR)];
fig.YLim=[1e-4 1e0];
fig.XTick=0:5:30;
saveas(h, ['SERvsSNR_TS_uniform(Delta=' num2str(Delta) ').eps'], 'epsc');
