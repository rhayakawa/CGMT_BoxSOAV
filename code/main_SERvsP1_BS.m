%%--------------------------------------------------------------------------------------------
% comparison between naive quantizer and asymptotically optimal quantizer
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
arrR=[-1 1];
L=length(arrR);
SNR=15;

arrP1=0.025:0.025:0.975;
arrSERvsP_NV=zeros(size(arrP1));
arrSERvsP_AO=zeros(size(arrP1));
arrSERvsP_NV_naive=zeros(size(arrP1));
arrSERvsP_AO_naive=zeros(size(arrP1));
for P1Index=1:length(arrP1)
  p1=arrP1(P1Index);
  disp(['p1=' num2str(p1)]);
  arrP=[p1 1-p1];

  if p1<0.5
    arrQ2=-0.24:0.005:0.02;
  elseif p1>0.5
    arrQ2=-0.02:0.005:0.24;
  else
    arrQ2=-0.04:0.005:0.04;
  end
  arrSER_theo_AO=zeros(size(arrQ2));
  arrSER_theo_NV=zeros(size(arrQ2));
  arrKappa=zeros(size(arrQ2));
  for Q2Index=1:length(arrQ2)
    Q2=arrQ2(Q2Index);
    arrCoef=[1 1-Q2];
    arrQ=[-Inf Q2 Inf];
    [~,arrSER_theo_NV(Q2Index),alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
    % threshold of AOQ
    arrKappa(Q2Index)=softThr(alpha_opt^(2)/(2*Delta)*log(p1/(1-p1)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
    arrThr=[-Inf arrKappa(Q2Index) Inf];
    arrSER_theo_AO(Q2Index)=SER_theo(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ,arrThr);
    if Q2==0
      arrSERvsP_NV_naive(P1Index)=arrSER_theo_NV(Q2Index);
      arrSERvsP_AO_naive(P1Index)=arrSER_theo_AO(Q2Index);
    end
  end
  arrSERvsP_NV(P1Index)=min(arrSER_theo_NV);
  arrSERvsP_AO(P1Index)=min(arrSER_theo_AO);
end

%% plot result
close all;
arrPlotIndeces=[1 4:4:36 39];
figure;
h=semilogy(arrP1,arrSERvsP_NV_naive,'--xk','LineWidth',1,'MarkerSize',10,'MarkerIndices',arrPlotIndeces);
hold on;
h=semilogy(arrP1,arrSERvsP_NV,'--ok','LineWidth',1,'MarkerSize',10,'MarkerIndices',arrPlotIndeces);
h=semilogy(arrP1,arrSERvsP_AO_naive,'-xk','LineWidth',1,'MarkerSize',10,'MarkerIndices',arrPlotIndeces);
h=semilogy(arrP1,arrSERvsP_AO,'-ok','LineWidth',1,'MarkerSize',10,'MarkerIndices',arrPlotIndeces);
grid on;

setLegend={'$\mathcal{Q}_{\mathrm{NV}} (\cdot)$ $(Q_{2}=0)$','$\mathcal{Q}_{\mathrm{NV}} (\cdot)$ (optimal)','$\mathcal{Q}_{\mathrm{AO}} (\cdot)$ $(Q_{2}=0)$','$\mathcal{Q}_{\mathrm{AO}} (\cdot)$ (optimal)'};
objLegend=legend(setLegend,'Location','best');
objLegend.Interpreter='latex';
objLegend.FontSize=16;
legend('boxoff');

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('$p_{1}$');
ylabel('asymptotic SER');
fig.XLim=[0 1];
fig.YLim=[1e-6 1e-2];
saveas(h, ['SERvsP1(Delta=' num2str(Delta) ',SNR=' num2str(SNR) ').eps'], 'epsc');
