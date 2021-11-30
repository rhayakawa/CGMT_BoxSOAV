%%--------------------------------------------------------------------------------------------
% comparison between naive quantizer and asymptotically optimal quantizer
% in discrete-valued sparse vector reconstruction
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
p0=0.9;
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
L=length(arrR);
SNR=15;

arrDelta=0.3:0.01:0.6;
arrSERvsDelta_nearest=zeros(size(arrDelta));
arrSERvsDelta_opt=zeros(size(arrDelta));
for DeltaIndex=1:length(arrDelta)
  Delta=arrDelta(DeltaIndex);
  disp(['Delta=' num2str(Delta)]);

  arrq2=[0.0001 0.001 0.002:0.002:0.04];
  arrSER_theo_opt=zeros(size(arrq2));
  arrSER_theo_nearest=zeros(size(arrq2));
  for q2Index=1:length(arrq2)
    q2=arrq2(q2Index);
    arrCoef=[1 q2 1];
    arrQ=arrCoef*[-1  1  1  1;...
                  -1 -1  1  1;...
                  -1 -1 -1  1;];
    arrQ(1)=-Inf;
    arrQ(L+1)=Inf;
    [~,arrSER_theo_nearest(q2Index),alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
    % threshold of AOQ
    kappa2=softThr(-1/2+alpha_opt^(2)/Delta*log((1-p0)/(2*p0)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
    arrThr=[-Inf kappa2 -kappa2 Inf];
    arrSER_theo_opt(q2Index)=SER_theo(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ,arrThr);
  end
  [arrSERvsDelta_nearest(DeltaIndex),tmp1]=min(arrSER_theo_nearest);
  [arrSERvsDelta_opt(DeltaIndex),tmp2]=min(arrSER_theo_opt);
end

%% plot result
close all;
figure;
h=semilogy(arrDelta,arrSERvsDelta_nearest,'--k','LineWidth',1);
hold on;
h=semilogy(arrDelta,arrSERvsDelta_opt,'-k','LineWidth',1);
grid on;

setLegend={'$\mathcal{Q}_{\mathrm{NV}} (\cdot)$','$\mathcal{Q}_{\mathrm{AO}} (\cdot)$'};
objLegend=legend(setLegend,'Location','northeast');
objLegend.Interpreter='latex';
objLegend.FontSize=16;
legend('boxoff');

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('$\Delta$');
ylabel('asymptotic SER');
fig.XLim=[min(arrDelta) max(arrDelta)];
fig.YLim=[1e-6 1e-1];
saveas(h, ['SERvsDelta_TS_NVAO(p0=' num2str(p0) ',SNR=' num2str(SNR) ').eps'], 'epsc');
