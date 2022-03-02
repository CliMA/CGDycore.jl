function cF=AverageFB(c)
cF=zeros(size(c,1),size(c,2),size(c,4),size(c,5)+1);
cF(:,:,:,1)=c(:,:,1,:,1);
cF(:,:,:,2:end-1)=0.5*(c(:,:,end,:,1:end-1)+c(:,:,1,:,2:end));
cF(:,:,:,end)=c(:,:,2,:,end);
end