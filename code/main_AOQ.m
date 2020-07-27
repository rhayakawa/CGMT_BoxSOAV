%%--------------------------------------------------------------------------------------------
% asymptotic probability density in binary vector reconstruction
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
Delta=0.6;
p1=0.3;
arrP=[p1 1-p1];
arrR=[-1 1];
L=length(arrR);
arrCoef=[0.5 0.5];
arrQ=arrCoef*[-1  1  1 ;...
              -1 -1  1 ;];
SNR=15;

%% theoretical asymptotic distribution
[~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
arrX=min(arrR):0.001:max(arrR);
for l=1:L
  arrX(arrX==arrR(l))=Inf;
end
arrX=arrX(arrX<Inf);
matp_x_hat_theo=zeros(L,length(arrX));
for l=1:L
  matp_x_hat_theo(l,:)=sqrt(Delta)/alpha_opt*arrP(l)*normpdf(sqrt(Delta)/alpha_opt*(softThr_inv(arrX,alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR)-arrR(l)));
end
arrP_X_hat_theo=zeros(1,L);
arrQ(1)=-Inf;
arrQ(L+1)=Inf;
for l=1:L
  for k=1:L
    arrP_X_hat_theo(l)=arrP_X_hat_theo(l)+arrP(k)*(normcdf(sqrt(Delta)/alpha_opt*(arrR(l)-arrR(k))+arrQ(l+1)/beta_opt)-normcdf(sqrt(Delta)/alpha_opt*(arrR(l)-arrR(k))+arrQ(l)/beta_opt));
  end
end

% threshold of AOQ
kappa=softThr(alpha_opt^(2)/(2*Delta)*log(p1/(1-p1)),alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR);

%% plot histogram
close all;
figure;
for l=1:L
  h=plot(arrX,matp_x_hat_theo(l,:),':k','LineWidth',1);
  hold on;
end
h=plot(arrX,sum(matp_x_hat_theo),'-k','LineWidth',1);
h=stem(arrR,arrP_X_hat_theo,'-^k','LineWidth',1);
h=plot([kappa kappa],[0 1],'--k','LineWidth',1);
grid on;

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('$\hat{x}$');
ylabel('probability density');
fig.XLim=[min(arrR)-0.1 max(arrR)+0.1];
fig.YLim=[0 0.6];
annotation('textarrow',[0.345,0.365],[0.45,0.305],'String','$p_{1} p_{\hat{X} \mid X=r_{1}} (\hat{x}) $','FontSize',20,'Interpreter','latex');
annotation('textarrow',[0.62,0.545],[0.27,0.3],'String','$p_{2} p_{\hat{X} \mid X=r_{2}} (\hat{x}) $','FontSize',20,'Interpreter','latex');
annotation('textbox',[0.45 0.03 0.1 0.1],'String','$\kappa_{2}^{\ast}$','FontSize',20,'LineStyle','none','Interpreter','latex');
saveas(h, ['AOQ(Delta=' num2str(Delta) ',p1=' num2str(p1) ',q1=' num2str(arrCoef(1)) ',q2=' num2str(arrCoef(2)) ',SNR=' num2str(SNR) ').eps'], 'epsc');
