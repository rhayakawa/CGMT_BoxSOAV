function value=Ex_f_prox(alpha,beta,Delta,arrP,arrR,arrCoef,arrQ)
  
  L=length(arrR);
  R_cap=zeros(L,L,L+1);
  for l=1:L
    R_cap(l,:,:)=R_cap(l,:,:)-sqrt(Delta)/alpha*arrR(l);
  end
  for l_prime=1:L
    R_cap(:,l_prime,:)=R_cap(:,l_prime,:)+sqrt(Delta)/alpha*arrR(l_prime);
  end
  for k=1:L+1
    R_cap(:,:,k)=R_cap(:,:,k)+1/beta*arrQ(k);
  end
  phi=normpdf(R_cap);
  Phi=normcdf(R_cap);
  
  arrIntegral=zeros(1,L);
  for l=1:L
    tmp=sum(arrCoef.*abs(arrR(1)-arrR));
    arrIntegral(l)=tmp*Phi(l,1,2);
    for k=2:(L-1)
      tmp=sum(arrCoef.*abs(arrR(k)-arrR));
      arrIntegral(l)=arrIntegral(l)+tmp*(Phi(l,k,k+1)-Phi(l,k,k));
    end
    tmp=sum(arrCoef.*abs(arrR(L)-arrR));
    arrIntegral(l)=arrIntegral(l)+tmp*(1-Phi(l,L,L));
    for k=2:L
      arrR_cap_tmp=R_cap(l,:,k);
      arrCoefR_cap=arrCoef.*arrR_cap_tmp;
      arrIntegral(l)=arrIntegral(l)+alpha/sqrt(Delta)*(arrQ(k)*(-phi(l,k,k)+phi(l,k-1,k))+(-sum(arrCoefR_cap(1:(k-1)))+sum(arrCoefR_cap(k:L)))*(Phi(l,k,k)-Phi(l,k-1,k)));
    end
  end
  
  value=arrIntegral*arrP';
  
end
