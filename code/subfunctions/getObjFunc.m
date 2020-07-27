function objFunc=getObjFunc(alpha,beta,Delta,arrP,arrR,arrCoef,arrQ,sigma2_v)
  % getObjFunc: compute the objective function of scalar optimization problem
  
  objFunc=alpha*beta*sqrt(Delta)/2+sigma2_v*beta*sqrt(Delta)/(2*alpha)-beta^(2)/2-alpha*beta/(2*sqrt(Delta));
  Ex1=Ex_f_prox(alpha,beta,Delta,arrP,arrR,arrCoef,arrQ);
  Ex2=Ex_prox2(alpha,beta,Delta,arrP,arrR,arrQ);
  objFunc=objFunc+Ex1+beta*sqrt(Delta)/(2*alpha)*Ex2;
  
end
