function eta=softThr_noBox(u,gamma,arrQ,arrR)
  % softThr_noBox: soft thresholding function for SOAV optimization without box constraint

  arrFoo=gamma*arrQ;

  L=length(arrR);
  eta=u-arrFoo(1);
  for l=1:L
    index=(u>arrR(l)+arrFoo(l));
    eta(index)=arrR(l);
    index2=(u>arrR(l)+arrFoo(l+1));
    eta(index2)=u(index2)-arrFoo(l+1);
  end
end

