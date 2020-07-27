function eta=softThr(u,gamma,arrQ,arrR)
  % softThr: soft thresholding function for Box-SOAV optimization

  arrFoo=gamma*arrQ;

  L=length(arrR);
  eta=arrR(1)*ones(size(u));
  for l=2:L
    index=(u>arrR(l-1)+arrFoo(l));
    eta(index)=u(index)-arrFoo(l);
    index2=(u>arrR(l)+arrFoo(l));
    eta(index2)=arrR(l);
  end
end

