function Average!(cA,c)
cA=reshape(.5*(c[:,:,2,:,:]+c[:,:,1,:,:]),
  size(c,1),size(c,2),size(c,4),size(c,5));
end
