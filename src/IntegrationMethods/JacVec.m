function k=JacVec(J,v)
n1=size(v,1);
n2=size(v,2);
n=n1*n2;
rRho=reshape(permute(v(:,:,1),[2,1]),n,1);
rTh=reshape(permute(v(:,:,5),[2,1]),n,1);
rw=reshape(permute(v(:,:,4),[2,1]),n,1);
sRho=J.JRhoW*rw;
sTh=J.JThW*rw;
sw=J.JWRho*rRho+J.JWTh*rTh;
k=zeros(size(v));
k(:,:,1)=permute(reshape(sRho,n2,n1),[2,1]);
k(:,:,4)=permute(reshape(sw,n2,n1),[2,1]);
k(:,:,5)=permute(reshape(sTh,n2,n1),[2,1]);
end