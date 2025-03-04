
function F=FcnNonLin(VSp,Grid,Param)
% VSp --> VCart
VCart = VSp2VCart(VSp,Param.rot);
FB=RiemannNonLin(VCart,Grid,Param);
FV=FluxVolumeNonLin(VCart,Grid,Param);
FCart=ResDG(FV,FB,Grid,Param);
%FCart --> FSp
FSp=VCart2VSp(FCart,Param.rot);
F=FSp+Source(VSp,Param);
end

function VCart=VSp2VCart(VSp,rot)
NX=size(VSp,2);
NV=size(VSp,1);
NF=size(VSp,4);
VCart=zeros(NV+1,NX,NX,NF);
VCart(1,:,:,:)=VSp(1,:,:,:);
for iF=1:NF
  VCart(2,:,:,iF)=rot(:,:,1,1,iF).*reshape(VSp(2,:,:,iF),NX,NX)...
    +rot(:,:,2,1,iF).*reshape(VSp(3,:,:,iF),NX,NX);
  VCart(3,:,:,iF)=rot(:,:,1,2,iF).*reshape(VSp(2,:,:,iF),NX,NX)...
    +rot(:,:,2,2,iF).*reshape(VSp(3,:,:,iF),NX,NX);
  VCart(4,:,:,iF)=rot(:,:,1,3,iF).*reshape(VSp(2,:,:,iF),NX,NX)...
    +rot(:,:,2,3,iF).*reshape(VSp(3,:,:,iF),NX,NX);
end
end
function VSp=VCart2VSp(VCart,rot)
NX=size(VCart,2);
NV=size(VCart,1);
NF=size(VCart,4);
VSp=zeros(NV-1,NX,NX,NF);
VSp(1,:,:,:)=VCart(1,:,:,:);
for iF=1:NF
  VSp(2,:,:,iF)=rot(:,:,1,1,iF).*reshape(VCart(2,:,:,iF),NX,NX)...
    +rot(:,:,1,2,iF).*reshape(VCart(3,:,:,iF),NX,NX)...
    +rot(:,:,1,3,iF).*reshape(VCart(4,:,:,iF),NX,NX);
  VSp(3,:,:,iF)=rot(:,:,2,1,iF).*reshape(VCart(2,:,:,iF),NX,NX)...
    +rot(:,:,2,2,iF).*reshape(VCart(3,:,:,iF),NX,NX)...
    +rot(:,:,2,3,iF).*reshape(VCart(4,:,:,iF),NX,NX);
end
end
