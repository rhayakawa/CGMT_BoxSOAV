function SER=SER_theo_nearest(alpha,beta,Delta,arrP,arrR,arrQ)
  % SER_theo_nearest: compute asymptotic symbol error rate with naive quantizer

  L=length(arrR);
  
  SuccessRate=0;
  SuccessRate=SuccessRate+arrP(1)*normcdf(sqrt(Delta)/alpha*( (-arrR(1)+arrR(2))/2 + alpha/(beta*sqrt(Delta))*arrQ(2) ));
  for l=2:(L-1)
    int1=sqrt(Delta)/alpha*( (arrR(l-1)-arrR(l))/2 + alpha/(beta*sqrt(Delta))*arrQ(l) );
    int2=sqrt(Delta)/alpha*( (-arrR(l)+arrR(l+1))/2 + alpha/(beta*sqrt(Delta))*arrQ(l+1) );
    SuccessRate=SuccessRate+arrP(l)*(normcdf(int2)-normcdf(int1));
  end
  SuccessRate=SuccessRate+arrP(L)*(1-normcdf(sqrt(Delta)/alpha*( (arrR(L-1)-arrR(L))/2 + alpha/(beta*sqrt(Delta))*arrQ(L) )));
  
  SER=1-SuccessRate;
  
end
