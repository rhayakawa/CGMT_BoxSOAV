function value=MSE_theo(alpha,beta,Delta,arrP,arrR,arrQ)
  % MSE_theo: compute asymptotic mean square error
  
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
    arrIntegral(l)=(arrR(1)-arrR(l))^(2)*Phi(l,1,2);
    for k=2:(L-1)
      arrIntegral(l)=arrIntegral(l)+(arrR(k)-arrR(l))^(2)*(Phi(l,k,k+1)-Phi(l,k,k));
    end
    arrIntegral(l)=arrIntegral(l)+(arrR(L)-arrR(l))^(2)*(1-Phi(l,L,L));
    for k=2:L
      arrIntegral(l)=arrIntegral(l)+alpha^(2)/Delta*(-R_cap(l,k,k)*phi(l,k,k)+Phi(l,k,k)+2/beta*arrQ(k)*phi(l,k,k)+1/beta^(2)*arrQ(k)^(2)*Phi(l,k,k));
    end
    for k=2:L
      arrIntegral(l)=arrIntegral(l)+alpha^(2)/Delta*(R_cap(l,k-1,k)*phi(l,k-1,k)-Phi(l,k-1,k)-2/beta*arrQ(k)*phi(l,k-1,k)-1/beta^(2)*arrQ(k)^(2)*Phi(l,k-1,k));
    end
  end
  
  value=arrIntegral*arrP';
  
end
