%%--------------------------------------------------------------------------------------------
% comparison between empirical performance and theoretical prediction 
% in binary sparse vector reconstruction
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
arrN=[50 100:100:1000];
Delta=0.7;
p1=0.8;
arrP=[p1 1-p1];
arrR=[0 1];
SNR=15;
nIteration=200;
nSample=100;

% parameter selection (asymptotically optimal)
arrQ2=-0.02:0.005:0.4;
arrSER_theo_AO=zeros(size(arrQ2));
arrKappa=zeros(size(arrQ2));
for Q2Index=1:length(arrQ2)
  Q2=arrQ2(Q2Index);
  arrCoef=[1 1-Q2];
  arrQ=[-Inf Q2 Inf];
  [~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
  % threshold of AOQ
  arrKappa(Q2Index)=softThr((arrR(1)+arrR(2))/2+alpha_opt^(2)/((arrR(2)-arrR(1))*Delta)*log(p1/(1-p1)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);
  arrThr=[-Inf arrKappa(Q2Index) Inf];
  arrSER_theo_AO(Q2Index)=SER_theo(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ,arrThr);
end
[minSER_theo,optQ2Index]=min(arrSER_theo_AO);
optQ2=arrQ2(optQ2Index);
optKappa=arrKappa(optQ2Index);

arrCoef=[1 1-optQ2];
arrQ=arrCoef*[-1 1 1; -1 -1 1;];
arrThr=[-Inf optKappa Inf];

arrSERvsN=zeros(1,length(arrN));
for NIndex=1:length(arrN)
  N=arrN(NIndex);
  disp(['N=' num2str(N)]);
  
  [MSE,SER]=get_empirical(N,Delta,arrP,arrR,arrQ,arrThr,SNR,nIteration,nSample);

  arrSERvsN(NIndex)=SER;
end

%% plot result
close all;
figure;
h=semilogy(arrN,arrSERvsN,'dk','LineWidth',1,'MarkerSize',10);
set(h, 'MarkerFaceColor', get(h,'Color'));
hold on;
h=semilogy([0 max(arrN)],[minSER_theo minSER_theo],'-k','LineWidth',1);
grid on;

xlabel('$N$');
ylabel('SER');
arrLegend={'empirical','theoretical'};
objLegend=legend(arrLegend);
objLegend.Interpreter='latex';
objLegend.Location='northeast';
objLegend.FontSize=20;
legend('boxoff');
fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
axis([0 1000 1e-5 1]);

saveas(h, ['SERvsN_BS(Delta=' num2str(Delta) ',p1=' num2str(p1) ',SNR=' num2str(SNR) ').eps'], 'epsc');

