function F=ResDG(FV,FB,Grid,Param)
NX=size(FV,2);
nV=size(FV,1);
NF=size(FV,4);
F=FV;
for iF=1:NF

  iE=Grid.Faces(iF).E(1);
  F(:,:,1,iF)=F(:,:,1,iF)-FB(:,:,iE)/Param.MX(1);
  
  iE=Grid.Faces(iF).E(2);
  F(:,NX,:,iF)=reshape(reshape(F(:,NX,:,iF),nV,NX)...
    -FB(:,:,iE)/Param.MX(NX),nV,1,NX,1);
  
  iE=Grid.Faces(iF).E(3);
  F(:,:,NX,iF)=F(:,:,NX,iF)+FB(:,:,iE)/Param.MX(NX);
  
  iE=Grid.Faces(iF).E(4);
  F(:,1,:,iF)=reshape(reshape(F(:,1,:,iF),nV,NX)...
    +FB(:,:,iE)/Param.MX(1),nV,1,NX,1);
  for iv=1:nV
    F(iv,:,:,iF)=reshape(F(iv,:,:,iF),NX,NX)./Param.J(:,:,iF);
  end
end
end