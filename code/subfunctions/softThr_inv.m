function eta_inv=softThr_inv(u,gamma,arrQ,arrR)
  % softThr_inv: inverse of soft thresholding function for Box-SOAV optimziation

  arrFoo=gamma*arrQ;

  L=length(arrR);
  eta_inv=zeros(size(u));
  for l=1:(L-1)
    index=(u>arrR(l));
    eta_inv(index)=u(index)+arrFoo(l+1);
  end
end

