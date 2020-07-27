function value=regularizer(s,arrCoef,arrR)

  value=0;
  L=length(arrR);
  for index=1:L
    r=arrR(index);
    value=value+arrCoef(index)*abs(s-r);
  end
end

