function k=SchurSolve(v,J,fac,Param)
n1=size(v,1);
n2=size(v,2);
n=n1*n2;
rRho=reshape(permute(v(:,:,1),[2,1]),n,1);
rTh=reshape(permute(v(:,:,5),[2,1]),n,1);
rw=reshape(permute(v(:,:,4),[2,1]),n,1);
invfac=1/fac;
invfac2=invfac/fac;
if Param.Damping
  sw=(spdiags(repmat(invfac2,n,1),0,n,n)-invfac*J.JWW-J.JWRho*J.JRhoW-J.JWTh*J.JThW)...
    \(invfac*rw+J.JWRho*rRho+J.JWTh*rTh);
else
  sw=(spdiags(repmat(invfac2,n,1),0,n,n)-J.JWRho*J.JRhoW-J.JWTh*J.JThW)...
    \(invfac*rw+J.JWRho*rRho+J.JWTh*rTh);
end
sRho=fac*(rRho+J.JRhoW*sw);
sTh=fac*(rTh+J.JThW*sw);
k=zeros(size(v));
k(:,:,1)=permute(reshape(sRho,n2,n1),[2,1]);
k(:,:,2:3)=fac*v(:,:,2:3);
k(:,:,4)=permute(reshape(sw,n2,n1),[2,1]);
k(:,:,5)=permute(reshape(sTh,n2,n1),[2,1]);
end