function V=RungeKuttaHEVI(V,fSlow,dt,FcnSlow,FcnFast,L,U,p,q,Param)
Vn=V;
RK=Param.RK;
NN=Param.NN;
nV=Param.nV;
fVE=zeros([size(V) RK.nStage]);
fVI=zeros([size(V) RK.nStage]);
for iStage=1:RK.nStage
  V=Vn;
  if iStage>1
    for jStage=1:iStage-1
      V=V+dt*RK.ARKE(iStage,jStage)*fVE(:,:,:,:,jStage)+...
          dt*RK.ARKI(iStage,jStage)*fVI(:,:,:,:,jStage);
    end
    v=reshape(V,NN*nV,1);
    v(q)=U\(L\v(p));
    V=reshape(v,size(Vn));
  end
  fVE(:,:,:,:,iStage)=FcnSlow(V,Param)+fSlow;
  fVI(:,:,:,:,iStage)=FcnFast(V,Param);
end
V=Vn;
for iStage=1:RK.nStage
  V=V+dt*RK.bRKE(iStage)*fVE(:,:,:,:,iStage)+...
      dt*RK.bRKI(iStage)*fVI(:,:,:,:,iStage);
end
end

