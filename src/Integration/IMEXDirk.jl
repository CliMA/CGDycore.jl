mutable struct CacheIMEXDirkStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  fV::AT4
  R::AT4
  dZ::AT4
  Z::AT5
  fEY::AT5
end

function Cache(backend,FT,IntMethod::IMEXDirkMethod,FE,M,nz,NumV)
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  R = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  dZ = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  Z = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  fEY = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  return CacheIMEXDirkStruct{FT,
                        typeof(Vn),
                        typeof(Z)}(
    Vn,
    fV,
    R,
    dZ,
    Z,
    fEY,
  )
end
#function IMEXSchur!(V,dt,FcnE,FcnI,Jac,CG,Global,Param)
function TimeIntegration!(IMEX::IMEXDirkMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  nStage=IMEX.nStage;
  fEY = Cache.fEY
  Z = Cache.Z
  fV = Cache.fV
  gammaD = IMEX.AI[2,2]
  R = Cache.R
  dZ = Cache.dZ
  Vn = Cache.Vn
  @views VnI = Vn[:,:,1:size(fV,3),1:size(fV,4)]
  dtau, = dt
  FcnE, FcnI = Fcn

# @views FcnE(fEY[:,:,:,1],V,CG,Global,Param)
  @. VnI = V
  @views FcnE(fEY[:,:,:,:,1],Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
  Jac(V,dtau*gammaD,FE,Metric,Phys,Aux,JCache,Global,DiscType)
  @views @. Z[:,:,:,:,1] = 0.0
  @inbounds for iStage = 2 : nStage 
    @. R = IMEX.D[iStage,1] * dtau * fEY[:,:,:,:,1]
    @inbounds for jStage = 2 : iStage - 1
      @views @. R = R + IMEX.D[iStage,jStage] * dtau * fEY[:,:,:,:,jStage] + 
        IMEX.E[iStage,jStage] * Z[:,:,:,:,jStage]
    end  
    @views @. Z[:,:,:,:,iStage] = Z[:,:,:,:,iStage-1]
    # Newton Iteration
    fac = gammaD * dtau
    @inbounds for iTer = 1 : 2
      @views @. VnI = V + Z[:,:,:,:,iStage]
      FcnI(fV,Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
      @views @. fV = (1 / fac) * (Z[:,:,:,:,iStage] - R) - fV
      Solve!(dZ,fV,JCache,dtau*gammaD,FE,Metric,Global,DiscType)
      @views @. Z[:,:,:,:,iStage] = Z[:,:,:,:,iStage] - dZ 
    end  
    @views @. VnI = V + Z[:,:,:,:,iStage]
    @views FcnE(fEY[:,:,:,:,iStage],Vn,FE,Metric,Phys,Aux,Exchange,Global,DiscType)
  end  
  @views @. V = V + IMEX.d[1] * dtau *  fEY[:,:,:,:,1]  
  @inbounds for jStage = 2 : nStage
    @views @. V = V + IMEX.d[jStage] * dtau *  fEY[:,:,:,:,jStage] + 
      IMEX.e[jStage] * Z[:,:,:,:,jStage]  
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

