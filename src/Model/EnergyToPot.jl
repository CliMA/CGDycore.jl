function Th=EnergyToPot(U,CG,Param)
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
W=[zeros(size(U,1),1)...
   ,U(:,:,wPos)];
WC=0.5*(W(:,1:nz)+W(:,2:nz+1));
KE=0.5*(U(:,:,uPos).*U(:,:,uPos)+U(:,:,vPos).*U(:,:,vPos)+WC.*WC);
Temp=(U(:,:,ThPos)-U(:,:,RhoPos).*(KE...
  +Param.Grav*repmat(Param.Grid.zP,1,size(U,1))'))/Param.Cvd;
Th=Temp.*(Param.p0./(Param.Rd.*Temp)).^(Param.Rd/Param.Cpd);
end