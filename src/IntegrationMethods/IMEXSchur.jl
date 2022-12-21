function IMEXSchur!(V,dt,FcnE,FcnI,Jac,CG,Global,Param)
  IMEX=Global.IMEX;
  nV1=size(V,1);
  nV2=size(V,2);
  nV3=size(V,3);
  nJ=nV1*nV2*nV3;
  nStage=IMEX.nStage;
  fEY=Global.Cache.Y
  Z=Global.Cache.Z
  fV=Global.Cache.fV
  R=Global.Cache.R
  dZ=Global.Cache.dZ
  Vn=Global.Cache.Vn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr
  Global.J.CompTri = true
  Global.J.CompJac = true
  J = Global.J

  @. Vn = V
  @views FcnE(fEY[:,:,:,1],V,CG,Global,Param)
  Jac(J,V,CG,Global,Param)
  @views @. Z[:,:,:,1] = 0.0
  @inbounds for iStage = 2 : nStage 
    @. R = 0.0
    @inbounds for jStage = 1 : iStage - 1
      @views @. R = R + IMEX.D[iStage,jStage] * dt * fEY[:,:,:,jStage] + 
        IMEX.E[iStage,jStage] * Z[:,:,:,jStage]
    end  
    @views @. Z[:,:,:,iStage] = Z[:,:,:,iStage-1]
    # Newton Iteration
    fac = IMEX.AI[iStage,iStage] * dt
    @inbounds for iTer = 1 : 2
      FcnI(fV,V,CG,Global,Param)  
      @views @. fV = Z[:,:,:,iStage] - R - fac  * fV
      SchurSolve!(dZ,fV,J,fac,Global)
      @views @. Z[:,:,:,iStage] = Z[:,:,:,iStage] - (1.0/fac) * dZ 
      @views @. V = Vn + Z[:,:,:,iStage]
    end  
    @views FcnE(fEY[:,:,:,iStage],V,CG,Global,Param)
  end  
  @. V = Vn
  @inbounds for jStage = 1 : nStage
    @views @. V = V + IMEX.d[jStage] * dt *  fEY[:,:,:,jStage] + 
      IMEX.e[jStage] * Z[:,:,:,jStage]  
  end
end  


function IMEXSchurCA!(V,dt,FcnE,FcnI,Jac,CG,Global,Param)
  IMEX=Global.IMEX;
  nV1=size(V,1);
  nV2=size(V,2);
  nV3=size(V,3);
  nJ=nV1*nV2*nV3;
  nStage=IMEX.nStage;
  fEY=Global.Cache.Y
  fIY=Global.Cache.Z
  fV=Global.Cache.fV
  R=Global.Cache.R
  dZ=Global.Cache.dZ
  Vn=Global.Cache.Vn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr
  Global.J.CompTri = true
  Global.J.CompJac = true
  J = Global.J

  @. Vn = V
  @views FcnE(fEY[:,:,:,1],V,CG,Global,Param)
  Jac(J,V,CG,Global,Param)
  @inbounds for iStage = 2 : nStage 
    @. V = Vn
    @inbounds for jStage = 1 : iStage - 1
      @views @. V = V + IMEX.AE[iStage,jStage] * dt * fEY[:,:,:,jStage] 
    end  
    @inbounds for jStage = 2 : iStage - 1
      @views @. V = V + IMEX.AI[iStage,jStage] * dt * fIY[:,:,:,jStage] 
    end  
    # Newton Iteration
    fac = IMEX.AI[iStage,iStage] * dt
    @. R = V
    @inbounds for iTer = 1 : 2
      FcnI(fV,V,CG,Global,Param)  
      @views @. fV = V - R - fac  * fV
      SchurSolve!(dZ,fV,J,fac,Global)
      @views @. V = V - (1.0/fac) * dZ 
    end  
    @. fIY[:,:,:,iStage] = (V - R) / fac
    @views FcnE(fEY[:,:,:,iStage],V,CG,Global,Param)
  end  
  @. V = Vn
  @inbounds for jStage = 2 : nStage
    @views @. V = V + IMEX.bE[jStage] * dt *  fEY[:,:,:,jStage] + 
      IMEX.bI[jStage] * dt * fIY[:,:,:,jStage]  
  end
end  

