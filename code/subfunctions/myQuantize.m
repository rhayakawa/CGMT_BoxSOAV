function value=myQuantize(x,arrR,arrThr)
  value=arrR(1)*ones(size(x));
  for index=2:length(arrR)
    value(x>arrThr(index))=arrR(index);
  end
end