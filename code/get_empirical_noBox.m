function [MSE,SER]=get_empirical_noBox(N,Delta,arrP,arrR,arrQ,arrThr,SNR,nIteration,nSample)
  % get_empirical_noBox: get empirical performance of SOAV optimization without box constraint
  %
  % Input
  %   N: dimension of unknown vector
  %   Delta: measurement ratio
  %   arrP: probability distribution of unknown vector
  %   arrR: candidate of unknown variables
  %   arrQ: array for Q_k
  %   arrThr: array for thresholds in final quantization
  %   SNR: singal-to-noise ratio
  %   nIteration: number of iterations
  %   nSample: number of samples
  %
  % Output
  %   MSE: mean square error
  %   SER: symbol error rate


  % add path
  addpath subfunctions;

  % number of measurements 
  M=round(N*Delta);
  % noise variance
  sigma2_v=arrP*(arrR.^(2)).'/(10^(SNR/10));

  % cumulative distribution
  L=length(arrR);
  matOne=ones(L,L);
  arrCDF=arrP*triu(matOne);

  rng('shuffle');

  SumMSE=0;
  SumSER=0;
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
      z=invMat*(z_til+gamma*x_MF);
      z_til=z_til+theta*(softThr_noBox(2*z-z_til,gamma,arrQ,arrR)-z);
    end
    SumMSE=SumMSE+norm(z-x)^(2)/N;
    SumSER=SumSER+nnz(myQuantize(z,arrR,arrThr)-x)/N;

  end
  MSE=SumMSE/nSample;
  SER=SumSER/nSample;

end