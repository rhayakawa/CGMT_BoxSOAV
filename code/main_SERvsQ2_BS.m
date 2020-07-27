%%--------------------------------------------------------------------------------------------
% asymptotic SER versus Q2 in binary sparse vector reconstruction
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
p0=0.8;
arrP=[p0 1-p0];
arrR=[0 1];
arrQ2=-0.02:0.002:0.1;
SNR=15;

arrMSE_theo=zeros(1,length(arrQ2));
arrSER_theo=zeros(1,length(arrQ2));
for Q2Index=1:length(arrQ2)
  Q2=arrQ2(Q2Index);
  disp(['Q2=' num2str(Q2)]);
  arrCoef=[1 1-Q2];
  arrQ=arrCoef*[-1 1 1; -1 -1 1;];

  [MSE,SER,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
  
  arrMSE_theo(Q2Index)=MSE;
  arrSER_theo(Q2Index)=SER;
end


%% plot SER
close all;
figure;
h=semilogy(arrQ2,arrSER_theo,'-k','LineWidth',1);
hold on;
grid on;

[minSER_theo,optIndex]=min(arrSER_theo);
optQ2=arrQ2(optIndex);
h=semilogy([optQ2 optQ2],[1e-5 1],'--k');
annotation('textarrow',[0.30,0.362],[0.19,0.14],'String','$Q_{2}^{\ast}$','FontSize',20,'Interpreter','latex');

xlabel('$Q_{2}$');
ylabel('asymptotic SER');
fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XTick=-0.02:0.02:0.1;
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
axis([min(arrQ2) max(arrQ2) 1e-5 1]);
saveas(h, ['SERvsQ2_BS(Delta=' num2str(Delta) ',p1=' num2str(p0) ',SNR=' num2str(SNR) ').eps'], 'epsc');
