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
  for iStage = 2 : nStage 
    @. V = 0.0
    for jStage = 1 : iStage - 1
      @views @. V = V + IMEX.D[iStage,jStage] * dt * fEY[:,:,:,jStage] + 
        IMEX.E[iStage,jStage] * Z[:,:,:,jStage]
    end  
    @. Z[:,:,:,iStage] = Z[:,:,:,iStage-1]
#   @. Z[:,:,:,iStage] = V
    @views StepImplicit!(Z[:,:,:,iStage],V,FcnI,Vn,J,IMEX.AI[iStage,iStage] * dt,CG,Global,Param)
    @views @. V = Vn + Z[:,:,:,iStage]
    @views FcnE(fEY[:,:,:,iStage],V,CG,Global,Param)
  end  
  @. V = Vn
  for jStage = 1 : nStage
     @views @. V = V + IMEX.d[jStage] * dt *  fEY[:,:,:,jStage] + 
      IMEX.e[jStage] * Z[:,:,:,jStage]  
  end
end  

function StepImplicit!(Z,R,FcnI,Vn,J,fac,CG,Global,Param)
  V = Global.Cache.fV
  @views dZ = Global.Cache.Z[:,:,:,Global.IMEX.nStage+1]
  @views fV = Global.Cache.Z[:,:,:,Global.IMEX.nStage+2]
  @. V = Vn + Z
  FcnI(fV,V,CG,Global,Param) 
  @. fV = Z - R - fac  * fV
  SchurSolve!(dZ,fV,J,fac,Global)
  @. Z = Z - (1.0/fac) * dZ 
  @. V = Vn + Z
  FcnI(fV,V,CG,Global,Param) 
  @. fV = Z - R - fac  * fV
  SchurSolve!(dZ,fV,J,fac,Global)
  @. Z = Z - (1.0/fac) * dZ 
end

