function F=Source(V,Param)
uPos=Param.uPos;
vPos=Param.vPos;
Omega=Param.Omega;
F=zeros(size(V));
NX=size(V,2);
NF=size(V,4);
switch Param.GridType
  case  'Radial'
    fac=reshape(2.0*Omega*Param.X(:,:,3,:)./sqrt(Param.X(:,:,1,:).^2....
    +Param.X(:,:,2,:).^2+Param.X(:,:,3,:).^2),NX,NX,NF);
    F(uPos,:,:,:)=fac.*reshape(V(vPos,:,:,:),NX,NX,NF);
    F(vPos,:,:,:)=-fac.*reshape(V(uPos,:,:,:),NX,NX,NF);  
end