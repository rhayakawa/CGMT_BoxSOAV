function SER=SER_theo(alpha,beta,Delta,arrP,arrR,arrQ,arrThr)
  % SER_theo: compute asymptotic symbol error rate with asymptotically optimal quantizer
  
  L=length(arrR);
  
  SuccessRate=0;
  for l=1:L
    int1=sqrt(Delta)/alpha*(arrThr(l)-arrR(l)) + arrQ(l)/beta;
    int2=sqrt(Delta)/alpha*(arrThr(l+1)-arrR(l)) + arrQ(l+1)/beta;
    SuccessRate=SuccessRate+arrP(l)*(normcdf(int2)-normcdf(int1));
  end
  SER=1-SuccessRate;
  
end
