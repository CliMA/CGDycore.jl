function V=RungeKuttaExplicit(V,dt,Fcn,CG,Param)
Vn=V;
RK=Param.RK;
fV=zeros([size(V) RK.nStage]);
for iStage=1:RK.nStage
  V=Vn;
  for jStage=1:iStage-1
    V=V+dt*RK.ARKE(iStage,jStage)*fV(:,:,:,jStage);
  end
  fV(:,:,:,iStage)=Fcn(V,CG,Param);
end
V=Vn;
for iStage=1:RK.nStage
  V=V+dt*RK.bRKE(iStage)*fV(:,:,:,iStage);
end
end
