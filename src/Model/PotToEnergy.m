function E=PotToEnergy(U,CG,Param)
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
W=[zeros(size(U,1),1)...
   ,U(:,:,wPos)];
WC=0.5*(W(:,1:nz)+W(:,2:nz+1));
KE=0.5*(U(:,:,Param.uPos).*U(:,:,Param.uPos)+U(:,:,Param.vPos).*U(:,:,Param.vPos)+WC.*WC);
T=TFromThetaRho(U(:,:,ThPos),U(:,:,RhoPos),Param);
E=U(:,:,RhoPos).*(Param.Cvd*T+KE+Param.Grav*repmat(Param.Grid.zP,1,size(U,1))');
end