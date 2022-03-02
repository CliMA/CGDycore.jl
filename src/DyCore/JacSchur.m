function [J]=JacSchur(U,CG,Param)
nz=Param.Grid.nz;
dz=Param.Grid.dz;
nCol=size(U,1);
nJ=nCol*nz;
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
% W=[zeros(size(U,1),1)...
%    ,U(:,:,wPos)];

RhoCol=permute(U(:,:,RhoPos),[2 1]);
ThCol=permute(U(:,:,ThPos),[2 1]);
if strcmp(Param.Thermo,'Energy')
  W=[zeros(size(U,1),1)...
    ,U(:,:,wPos)];
  WC=0.5*(W(:,1:nz)+W(:,2:nz+1));
  KE=0.5*(U(:,:,Param.uPos).*U(:,:,Param.uPos)+U(:,:,Param.vPos).*U(:,:,Param.vPos)+WC.*WC);
  Pres=Pressure(U(:,:,ThPos),U(:,:,RhoPos),KE,Param);
  PresCol=permute(Pres,[2 1]);
  PhiCol=Param.Grav*repmat(Param.Grid.zP,1,size(U,1));
else
  Pres=Pressure(U(:,:,ThPos),U(:,:,ThPos),U(:,:,ThPos),Param);
  PresCol=permute(Pres,[2 1]);
end


if strcmp(Param.ModelType,'Curl')
  D=zeros(nz,nCol);
  D(1:nz-1,:)=-0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:))/dz;
  J.JRhoW=spdiags([reshape(-D,nJ,1) reshape(D,nJ,1)],[-1 0]...
    ,nJ,nJ);
  if strcmp(Param.Thermo,'Energy')
    dPdTh=(Param.Rd/Param.Cvd)./RhoCol;
  else
    dPdTh=dPresdTh(ThCol,Param);
  end
  Dp=zeros(nz,nCol);
  Dm=zeros(nz,nCol);
  Dp(1:nz-1,:)=dPdTh(2:nz,:)...
    ./(0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:)))...
    /dz;
  Dm(1:nz-1,:)=dPdTh(1:nz-1,:)...
    ./(0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:)))...
    /dz;
  J.JWTh=spdiags([[reshape(Dm,nJ,1);0] [0;reshape(-Dp,nJ,1)]],[0 1]...
    ,nJ,nJ);
  
  if strcmp(Param.Thermo,'Energy')
    Dp=zeros(nz,nCol);
    Dm=zeros(nz,nCol);
    Dp(1:nz-1,:)=(Param.Rd/Param.Cvd)...
      ./(0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:)))...
      .*(PhiCol(2:nz,:))/dz;
    Dm(1:nz-1,:)=(Param.Rd/Param.Cvd)...
      ./(0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:)))...
      .*(PhiCol(1:nz-1,:))/dz;
    D(1:nz-1,:)=-0.5*Param.Grav...
      ./(0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:)));
    J.JWRho=spdiags([[reshape(-Dm+D,nJ,1);0] [0;reshape(Dp+D,nJ,1)]],[0 1]...
      ,nJ,nJ);
  else
    D(1:nz-1,:)=0.5*(PresCol(2:nz,:)-PresCol(1:nz-1,:))/Param.Grid.dz...
      ./(0.5*(RhoCol(1:nz-1,:)+RhoCol(2:nz,:))).^2;
    J.JWRho=spdiags([[reshape(D,nJ,1);0] [0;reshape(D,nJ,1)]],[0 1]...
      ,nJ,nJ);
  end
  
  if strcmp(Param.Thermo,'Energy')
    D(1:nz-1,:)=-0.5*(ThCol(1:nz-1,:)+ThCol(2:nz,:)...
      +PresCol(1:nz-1,:)+PresCol(2:nz,:))/dz;
  else
    D(1:nz-1,:)=-0.5*(ThCol(1:nz-1,:)+ThCol(2:nz,:))/dz;
  end
  J.JThW=spdiags([reshape(-D,nJ,1) reshape(D,nJ,1)],[-1 0]...
    ,nJ,nJ);
  if Param.Damping
    K=permute(DampingKoeff(CG,Param),[2 1]);
    J.JWW=spdiags(reshape(K,nJ,1),0,nJ,nJ);
  end
else
  % Conservative
  D=zeros(nz,nCol);
  D(1:nz-1,:)=-1/dz;
  J.JRhoW=spdiags([reshape(-D,nJ,1) reshape(D,nJ,1)],[-1 0]...
    ,nJ,nJ);
  
  dPdTh=dPresdTh(ThCol,Param);
  Dp=zeros(nz,nCol);
  Dm=zeros(nz,nCol);
  Dp(1:nz-1,:)=dPdTh(2:nz,:)/dz;
  Dm(1:nz-1,:)=dPdTh(1:nz-1,:)/dz;
  J.JWTh=spdiags([[reshape(Dm,nJ,1);0] [0;reshape(-Dp,nJ,1)]],[0 1]...
    ,nJ,nJ);
  D(1:nz-1,:)=-0.5*Param.Grav;
  J.JWRho=spdiags([[reshape(D,nJ,1);0] [0;reshape(D,nJ,1)]],[0 1]...
    ,nJ,nJ);
  Temp=ThCol./RhoCol;
  D(1:nz-1,:)=-0.5*(Temp(1:nz-1,:)+Temp(2:nz,:))/dz;
  J.JThW=spdiags([reshape(-D,nJ,1) reshape(D,nJ,1)],[-1 0]...
    ,nJ,nJ);
end
end