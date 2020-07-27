function [MSE,SER,alpha_opt,beta_opt]=get_theoretical(Delta,arrP,arrR,arrCoef,arrQ,SNR)
  % get_theoretical: get optimal values for scalar max-min optimization
  %
  % Input
  %   Delta: measurement ratio
  %   arrP: probability distribution of unknown vector
  %   arrR: candidate of unknown variables
  %   arrCoef: coefficients in regularizer
  %   arrQ: array for Q_k
  %   SNR: singal-to-noise ratio
  %
  % Output
  %   MSE: asymptotic mean square error
  %   SER: asymptotic symbol error rate with naive quantizer
  %   alpha_opt: optimal value of alpha
  %   beta_opt: optimal value of beta


  % add path
  addpath('subfunctions');

  % noise variance
  sigma2_v=arrP*(arrR.^(2)).'/(10^(SNR/10));

  alphaIteration=100;
  betaIteration=100;

  % optimization of alpha and beta via ternary search
  alphaMin=0;
  alphaMax=10;
  for alphaIterationIndex=1:alphaIteration
    alpha1=(2*alphaMin+alphaMax)/3;
    alpha2=(alphaMin+2*alphaMax)/3;
    % objFunc for alpha1
    betaMin=0;
    betaMax=10;
    for betaIterationIndex=1:betaIteration
      beta1=(2*betaMin+betaMax)/3;
      beta2=(betaMin+2*betaMax)/3;
      objFunc1=getObjFunc(alpha1,beta1,Delta,arrP,arrR,arrCoef,arrQ,sigma2_v);
      objFunc2=getObjFunc(alpha1,beta2,Delta,arrP,arrR,arrCoef,arrQ,sigma2_v);
      if objFunc1>objFunc2
        betaMax=beta2;
      else
        betaMin=beta1;
      end
      if (beta2-beta1)<1e-6
        break;
      end
    end
    objFunc_alpha1=objFunc1;
    % objFunc for alpha2
    betaMin=0;
    betaMax=10;
    for betaIterationIndex=1:betaIteration
      beta1=(2*betaMin+betaMax)/3;
      beta2=(betaMin+2*betaMax)/3;
      objFunc1=getObjFunc(alpha2,beta1,Delta,arrP,arrR,arrCoef,arrQ,sigma2_v);
      objFunc2=getObjFunc(alpha2,beta2,Delta,arrP,arrR,arrCoef,arrQ,sigma2_v);
      if objFunc1>objFunc2
        betaMax=beta2;
      else
        betaMin=beta1;
      end
      if (beta2-beta1)<1e-6
        break;
      end
    end
    objFunc_alpha2=objFunc1;
    if objFunc_alpha1<objFunc_alpha2
      alphaMax=alpha2;
    else
      alphaMin=alpha1;
    end
    if (alpha2-alpha1)<1e-6
      break;
    end
  end
  alpha_opt=alpha1;
  beta_opt=beta1;

  MSE=MSE_theo(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ);
  SER=SER_theo_nearest(alpha_opt,beta_opt,Delta,arrP,arrR,arrQ);

end
