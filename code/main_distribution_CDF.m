%%--------------------------------------------------------------------------------------------
% comparison between empirical distribution and asymptotic distribution
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
Delta=0.75;
M=round(N*Delta);
p0=0.6;
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
q2=0.005;
arrCoef=[1 q2 1];
arrQ=arrCoef*[-1  1  1  1;...
              -1 -1  1  1;...
              -1 -1 -1  1;];
SNR=20;
nIteration=1000;

% noise variance
sigma2_v=arrP*(arrR.^(2)).'/(10^(SNR/10));

% cumulative distribution
L=length(arrR);
matOne=ones(L,L);
arrCDF=arrP*triu(matOne);

rng('shuffle');

%% empirical distribution
nSample=20;
X_hat=zeros(N,nSample);
for sampleIndex=1:nSample
  % unknown discrete-valued vector
  x_rand=rand(N,1);
  x=ones(N,1)*arrR(1);
  for valueIndex=2:L
    x(x_rand>=arrCDF(valueIndex-1))=arrR(valueIndex);
  end
  % measurement matrix
  A=randn(M,N)/sqrt(N);
  % additive noise vector
  v=randn(M,1)*sqrt(sigma2_v);
  % linear measurements
  y=A*x+v;

  gamma=1;
  invMat=(eye(N)+gamma*(A'*A))^(-1);
  x_MF=A'*y;

  % SOAV optimizaion via Douglas-Rachford algorithm
  theta=1.9;
  z=zeros(N,1);
  z_til=zeros(N,1);
  for k=2:nIteration
    z=softThr(z_til,gamma,arrQ,arrR);
    z_til=z_til+theta*(invMat*(2*z-z_til+gamma*x_MF)-z);
  end
  X_hat(:,sampleIndex)=z;
end

%% theoretical asymptotic distribution
[~,~,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR);
arrP_X_hat_theo=zeros(1,L);
arrQ(1)=-Inf;
arrQ(L+1)=Inf;
for l=1:L
  for k=1:L
    arrP_X_hat_theo(l)=arrP_X_hat_theo(l)+arrP(k)*(normcdf(sqrt(Delta)/alpha_opt*(arrR(l)-arrR(k))+arrQ(l+1)/beta_opt)-normcdf(sqrt(Delta)/alpha_opt*(arrR(l)-arrR(k))+arrQ(l)/beta_opt));
  end
end

arrX=min(arrR)+0.0001:0.001:max(arrR);
arrCDF_theo=zeros(size(arrX));
for l=1:L
  arrCDF_theo=arrCDF_theo+arrP(l)*normcdf(sqrt(Delta)/alpha_opt*(softThr_inv(arrX,alpha_opt/(beta_opt*sqrt(Delta)),arrQ,arrR)-arrR(l)));
end

%% plot histogram
close all;
figure;
arrEdges=[-10 -1:0.05:1 10 15];
h=histogram(X_hat,arrEdges,'LineStyle','-','EdgeColor','auto','DisplayStyle','stairs','Normalization','cdf','LineWidth',1);
h.EdgeColor='k';
hold on;
h=plot(arrX,arrCDF_theo,'--k','LineWidth',1);
grid on;

fig=gca;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
xlabel('$\hat{x}$');
ylabel('CDF');
fig.XLim=[min(arrR)-0.1 max(arrR)+0.1];
fig.YLim=[0 1];
saveas(h, ['distribution_CDF(N=' num2str(N) ',Delta=' num2str(Delta) ',p0=' num2str(p0) ',q2=' num2str(arrCoef(2)) ',SNR=' num2str(SNR) ').eps'], 'epsc');
