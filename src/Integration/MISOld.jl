function MultirateMIS!(V,dt,FcnE,FcnI,CG,Global,Param)
  MIS=Param.MIS
  OrdPoly=Param.OrdPolyX
  fVSlow=zeros([size(V) MIS.nStage])
  Th=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,MIS.nStage)
  Dp=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,MIS.nStage)
  pD=zeros(OrdPoly+1,OrdPoly+1,Param.Grid.NumFaces,MIS.nStage)
  ThE=zeros(OrdPoly+1,Param.Grid.NumEdges,MIS.nStage)
  vn=V
  VSlow=zeros([size(V) MIS.nStage+1])
  VSlow(:,:,:,:,1)=vn
  for iStage=1:MIS.nStage
    VSlow(:,:,:,:,iStage+1)=vn
    [fVSlow(:,:,:,:,iStage),Th(:,:,:,iStage),ThE(:,:,iStage),Dp(:,:,:,iStage),pD(:,:,:,iStage)]...
      =FcnSlowNonLin(VSlow(:,:,:,:,iStage),Param)
    fSlow=MIS.A(iStage+1,1)*fVSlow(:,:,:,:,iStage)
    Param.ThSlow=MIS.A(iStage+1,1)*Th(:,:,:,iStage)
    Param.DpSlow=MIS.A(iStage+1,1)*Dp(:,:,:,iStage)
    Param.pDSlow=MIS.A(iStage+1,1)*pD(:,:,:,iStage)
    Param.ThSlowE=MIS.A(iStage+1,1)*ThE(:,:,iStage)
    for jStage=2:iStage
      fSlow=fSlow+MIS.A(iStage+1,jStage)*fVSlow(:,:,:,:,jStage)
      fSlow=fSlow+MIS.G(iStage+1,jStage)*(VSlow(:,:,:,:,jStage)-vn)/dt
      Param.ThSlow=Param.ThSlow+MIS.A(iStage+1,jStage)*Th(:,:,:,jStage)
      Param.DpSlow=Param.DpSlow+MIS.A(iStage+1,jStage)*Dp(:,:,:,jStage)
      Param.pDSlow=Param.pDSlow+MIS.A(iStage+1,jStage)*pD(:,:,:,jStage)
      Param.ThSlowE=Param.ThSlowE+MIS.A(iStage+1,jStage)*ThE(:,:,jStage)
      VSlow(:,:,:,:,iStage+1)=VSlow(:,:,:,:,iStage+1)...
        +MIS.D(iStage+1,jStage)*(VSlow(:,:,:,:,jStage)-vn)
    end
    nSmall=ceil(MIS.d(iStage+1)*dt/Param.dtFast)
    dtau=MIS.d(iStage+1)*dt/nSmall
    for iSmall=1:nSmall
      switch  Param.MISFastType
        case 'RungeKuttaEx'
          VSlow(:,:,:,:,iStage+1)=RungeKutta(VSlow(:,:,:,:,iStage+1),...
             dtau,fSlow,Param)
        case 'RungeKuttaMIS'
          VSlow(:,:,:,:,iStage+1)=RungeKuttaMIS(VSlow(:,:,:,:,iStage+1),...
            dtau,fSlow,Param)  
        case 'RungeKuttaIMEX1'
          J=JacFcnFast1(Param)
          [L,U,p,q]=lu(speye(Param.NN*Param.nV)-dtImp*RK.gRKI*J,'vector')
          VSlow(:,:,:,:,iStage+1)=RungeKuttaImexV1(VSlow(:,:,:,:,iStage+1),...
            hSlow,hSlowE,dtau,fSlow,Param)
        case 'RungeKuttaIMEX2'
          VSlow(:,:,:,:,iStage+1)=RungeKuttaImexV2(VSlow(:,:,:,:,iStage+1),...
            hSlow,hSlowE,dtau,fSlow,Param)
      end
    end
  end
  V=VSlow(:,:,:,:,MIS.nStage+1)
end

function V=RungeKutta(V,dt,fSlow,Param)
  Vn=V
  RK=Param.RKFast
  fV=zeros([size(V) RK.nStage])
  for iStage=1:RK.nStage
    V=Vn
    for jStage=1:iStage-1
      V=V+dt*RK.ARKE(iStage,jStage)*fV(:,:,:,:,jStage)
    end
    fV(:,:,:,:,iStage)=FcnFastNonLin(V,Param)+fSlow
  end
  V=Vn
  for iStage=1:RK.nStage
    V=V+dt*RK.bRKE(iStage)*fV(:,:,:,:,iStage)
  end
end

function V=RungeKuttaImexV1(V,hSlow,hSlowE,dt,fSlow,Param)
  Vn=V
  RK=Param.RK
  NN=Param.NN
  nV=Param.nV
  fVE=zeros([size(V) RK.nStage])
  fVI=zeros([size(V) RK.nStage])
  for iStage=1:RK.nStage
    V=Vn
    if iStage>1
      for jStage=1:iStage-1
        V=V+dt*RK.ARKE(iStage,jStage)*fVE(:,:,:,:,jStage)+...
            dt*RK.ARKI(iStage,jStage)*fVI(:,:,:,:,jStage)
      end
      v=reshape(V,NN*nV,1)
      v(q)=U\(L\v(p))
      V=reshape(v,size(Vn))
    end
    fVE(:,:,:,:,iStage)=FcnFastNonLin1(V,hSlow,hSlowE,Param)+fSlow
    fVI(:,:,:,:,iStage)=FcnFastNonLin2(V,hSlow,hSlowE,Param)
  end
  V=Vn
  for iStage=1:RK.nStage
      V=V+dt*RK.bRKE(iStage)*fVE(:,:,:,:,iStage)+...
          dt*RK.bRKI(iStage)*fVI(:,:,:,:,iStage)
  end
end
function V=RungeKuttaImex2(V,hSlow,hSlowE,dt,fSlow,Param)
  Vn=V
  RK=Param.RK
  NN=Param.NN
  nV=Param.nV
  fVE=zeros([size(V) RK.nStage])
  fVI=zeros([size(V) RK.nStage])
  for iStage=1:RK.nStage
    V=Vn
    if iStage>1
      for jStage=1:iStage-1
        V=V+dt*RK.ARKE(iStage,jStage)*fVE(:,:,:,:,jStage)+...
            dt*RK.ARKI(iStage,jStage)*fVI(:,:,:,:,jStage)
      end
      v=reshape(V,NN*nV,1)
      v(q)=U\(L\v(p))
      V=reshape(v,size(Vn))
    end
    fVE(:,:,:,:,iStage)=FcnFastNonLin2(V,hSlow,hSlowE,Param)+fSlow
    fVI(:,:,:,:,iStage)=FcnFastNonLin1(V,hSlow,hSlowE,Param)
  end
  V=Vn
  for iStage=1:RK.nStage
    V=V+dt*RK.bRKE(iStage)*fVE(:,:,:,:,iStage)+...
        dt*RK.bRKI(iStage)*fVI(:,:,:,:,iStage)
  end
end



