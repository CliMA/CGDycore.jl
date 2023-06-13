function VerticalDiffusionScalarRho!(Fc,FRho,c,Rho,K,CG,dXdxI33,J,ThreadCache)
  @unpack TCacheCC1, TCacheCC2 = ThreadCache

  nz = size(Fc,3)
  gradqCG  = TCacheCC1[Threads.threadid()]
  qCG  = TCacheCC2[Threads.threadid()]

  @views @. qCG = c / Rho 
# Gradient computation
  @inbounds for iz = 1:nz-1
    @views @. gradqCG[:,:,iz+1] = (K[:,:,iz] + K[:,:,iz+1]) * (qCG[:,:,iz+1] - qCG[:,:,iz]) * 
      (dXdxI33[:,:,2,iz] + dXdxI33[:,:,1,iz+1]) / ( J[:,:,2,iz] + J[:,:,1,iz+1]) 
  end
# Divergence  
  @views @. Fc[:,:,1] += gradqCG[:,:,2] * dXdxI33[:,:,2,1]
  @views @. FRho[:,:,1] += gradqCG[:,:,2] * dXdxI33[:,:,2,1]
  @inbounds for iz = 2:nz-1
    @views @. Fc[:,:,iz] += (gradqCG[:,:,iz+1] * dXdxI33[:,:,2,iz] - gradqCG[:,:,iz] * dXdxI33[:,:,1,iz])
    @views @. FRho[:,:,iz] += (gradqCG[:,:,iz+1] * dXdxI33[:,:,2,iz] - gradqCG[:,:,iz] * dXdxI33[:,:,1,iz])
  end
  @views @. Fc[:,:,nz] -= gradqCG[:,:,nz] * dXdxI33[:,:,1,nz]
  @views @. FRho[:,:,nz] -= gradqCG[:,:,nz] * dXdxI33[:,:,1,nz]
end

function VerticalDiffusionScalar!(Fc,c,Rho,K,CG,dXdxI33,J,ThreadCache)
  @unpack TCacheCC1, TCacheCC2 = ThreadCache

  nz = size(Fc,3)
  gradqCG  = TCacheCC1[Threads.threadid()]
  qCG  = TCacheCC2[Threads.threadid()]

  @views @. qCG = c / Rho 
# Gradient computation
  @inbounds for iz = 1:nz-1
    @views @. gradqCG[:,:,iz+1] = (K[:,:,iz] + K[:,:,iz+1]) * (qCG[:,:,iz+1] - qCG[:,:,iz]) * 
      (dXdxI33[:,:,2,iz] + dXdxI33[:,:,1,iz+1]) / ( J[:,:,2,iz] + J[:,:,1,iz+1]) 
  end
# Divergence  
  @views @. Fc[:,:,1] += gradqCG[:,:,2] * dXdxI33[:,:,2,1]
  @inbounds for iz = 2:nz-1
    @views @. Fc[:,:,iz] += (gradqCG[:,:,iz+1] * dXdxI33[:,:,2,iz] - gradqCG[:,:,iz] * dXdxI33[:,:,1,iz])
  end
  @views @. Fc[:,:,nz] -= gradqCG[:,:,nz] * dXdxI33[:,:,1,nz]
end

function VerticalDiffusionMomentum!(FuC,FvC,uC,vC,Rho,K,dXdxI33,J,Cache)

  nz = size(uC,3)
  graduC = Cache.CacheC1
  gradvC = Cache.CacheC2
  @. K = 200.0  # noch zu berechnen
# Gradient computation
  @inbounds for iz = 1:nz-1
    @. graduC[:,:,iz+1] = 0.5 * (K[:,:,iz] + K[:,:,iz+1]) * (uC[:,:,iz+1] - uC[:,:,iz]) * 
      (dXdxI33[:,:,2,iz] + dXdxI33[:,:,1,iz+1]) / (J[:,:,2,iz] + J[:,:,1,iz+1])
    @. gradvC[:,:,iz+1] = 0.5 * (K[:,:,iz] + K[:,:,iz+1]) * (vC[:,:,iz+1] - vC[:,:,iz]) * 
      (dXdxI33[:,:,2,iz] + dXdxI33[:,:,1,iz+1]) / (J[:,:,2,iz] + J[:,:,1,iz+1])
  end
# Divergence
  @. FuC[:,:,1] += Rho[:,:,1] * graduC[:,:,2] * dXdxI33[:,:,2,1]
  @. FvC[:,:,1] += Rho[:,:,1] * gradvC[:,:,2] * dXdxI33[:,:,2,1]
  @inbounds for iz = 2:nz-1
    @. FuC[:,:,iz] += Rho[:,:,iz] * (graduC[:,:,iz+1] * dXdxI33[:,:,2,iz]  - graduC[:,:,iz] * dXdxI33[:,:,1,iz])
    @. FvC[:,:,iz] += Rho[:,:,iz] * (gradvC[:,:,iz+1] * dXdxI33[:,:,2,iz]  - gradvC[:,:,iz] * dXdxI33[:,:,1,iz])
  end
  @. FuC[:,:,nz] -= Rho[:,:,nz] * graduC[:,:,nz] * dXdxI33[:,:,1,nz]
  @. FvC[:,:,nz] -= Rho[:,:,nz] * gradvC[:,:,nz] * dXdxI33[:,:,1,nz]
end

function BoundaryFluxScalar!(Fc,c,cS,CG,Global,Param,iF)
  @. Fc -= Param.CTr * Global.Cache.uStar[:,:,iF] * (c - cS) *
    Global.Metric.dXdxIF[:,:,1,3,3,iF]
end

function BoundaryFluxMomentum!(FuC,FvC,uC,vC,w,Global,Param,iF)
 
  @views nS = Global.Metric.nS[:,:,:,iF]
  @views FS = Global.Metric.FS[:,:,iF]
  OP = size(uC,1)
  @inbounds for j = 1 : OP
    @inbounds for i = 1 : OP
      nSTV = uC[i,j] * nS[i,j,1] + vC[i,j] * nS[i,j,2] + w[i,j] * nS[i,j,3]
      uStar = sqrt(uC[i,j] * uC[i,j] + vC[i,j] * vC[i,j] + w[i,j] * w[i,j] - nSTV * nSTV) 
      FuC[i,j] -= Param.CMom * uStar * FS[i,j] * (uC[i,j] - nSTV * nS[i,j,1])
      FvC[i,j] -= Param.CMom * uStar * FS[i,j] * (vC[i,j] - nSTV * nS[i,j,2])
    end  
  end
end

function BoundaryFluxScalar!(Fc,Th,Rho,Tr,p,CG,Global,Param,iF)
  if Global.Model.Problem == "HeldSuarezMoistSphere" || Global.Model.Problem == "HeldSuarezMoistSphereOro"
    OP = CG.OrdPoly+1  
    ThPos=Global.Model.ThPos  
    RhoPos=Global.Model.RhoPos  
    RhoVPos=Global.Model.RhoVPos  
    NumV=Global.Model.NumV  
    @views TSurf = Global.Cache.TSurf[:,:,iF]
    @views uStar = Global.Cache.uStar[:,:,iF]
    CE = Param.CE
    CH = Param.CH
    Rd = Global.Phys.Rd
    Cpd = Global.Phys.Cpd
    Rv = Global.Phys.Rv
    Cpv = Global.Phys.Cpv
    p0 = Global.Phys.p0
    @views dXdxI = Global.Metric.dXdxI[:,:,1,1,3,3,iF]
    @inbounds for j = 1:OP
      @inbounds for i = 1:OP
       RhoV = max(Tr[i,j,RhoVPos], 0.0) 
       (FTh,FRho,FRhoV) = BoundaryFluxHeldSuarez(
         Th[i,j],Rho[i,j],RhoV,p[i,j],TSurf[i,j],dXdxI[i,j],CH,CE,uStar[i,j], Global.Phys)
       Fc[i,j,ThPos] += FTh
       Fc[i,j,RhoPos] += FRho
       Fc[i,j,RhoVPos+NumV] += FRhoV
      end   
    end   
  end    
end

function BoundaryFluxHeldSuarez(Th,Rho,RhoV,p,TSurf,dXdxIF,CH,CE,uStar, Phys)
  RhoD = Rho - RhoV
  Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
  Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV
  T = p / Rm
  p_vs = fpvs(TSurf,Phys)
  RhoVSurface = p_vs / (Phys.Rv * TSurf) 
  LatFlux = - 2.0 * CE * uStar * dXdxIF  * (RhoV - RhoVSurface) 
  SensFlux = - 2.0 * CH * uStar * dXdxIF  * (T - TSurf) 
  FRho = LatFlux
  FRhoV = LatFlux
  PrePi=(p / Phys.p0)^(Rm / Cpml)
  FTh = Th * (SensFlux / T + ((Phys.Rv / Rm) - 1.0 / Rho - log(PrePi)*(Phys.Rv / Rm - Phys.Cpv / Cpml)) *  LatFlux)  
  return (FTh,FRho,FRhoV)
end


function uStarCoefficient!(uStar,U,V,WC,CG,dXdxI,nS)
# Computation norm_v_a
# |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)
  OP = CG.OrdPoly+1  
  @inbounds for j = 1:OP
    @inbounds for i = 1:OP
      v1 = U[i,j] 
      v2 = V[i,j] 
      WS = -(dXdxI[i,j,3,1]* v1 + 
        dXdxI[i,j,3,2] * v2) / 
        dXdxI[i,j,3,3]
      w = 0.5 * (WS + WC[i,j])  
      uStar[i,j] = sqrt(v1 * v1 + v2 * v2 + w * w - 
        (nS[i,j,1] * v1 + nS[i,j,2] * v2 + nS[i,j,3] * w)^2)
    end
  end  
end

function Cd_coefficient!(CdTh,CdTr,CG,Global,Param,iF) 
  if Global.Model.Problem == "HeldSuarezMoistSphere" || Global.Model.Problem == "HeldSuarezMoistSphereOro"
    @views uStar = Global.Cache.uStar[:,:,iF]
    OP = CG.OrdPoly + 1  
    @inbounds for jP = 1 : OP
      @inbounds for iP = 1 : OP
        CdTh[iP,jP] = Param.CH * uStar[iP,jP]
        CdTr[iP,jP,1] = Param.CE * uStar[iP,jP] 
        CdTr[iP,jP,2] = 0
      end
    end  
  end    
end
function eddy_diffusivity_coefficient!(K,U,V,WC,Rho,p,CG,Global,Param,iF) 
  if Global.Model.Problem == "HeldSuarezMoistSphere" || Global.Model.Problem == "HeldSuarezMoistSphereOro"
    CE = Param.CE 
    p_pbl = Param.p_pbl 
    p_strato = Param.p_strato 
    OP = CG.OrdPoly + 1  
    nz = Global.Grid.nz
#   Computation norm_v_a  
#   |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)  
    @views uStar = Global.Cache.uStar[:,:,iF]
    @inbounds for jP = 1 : OP
      @inbounds for iP = 1 : OP
#       K[iP,jP,1] = CE * uStar[iP,jP] * (Global.Metric.J[iP,jP,1,1,iF] + Global.Metric.J[iP,jP,2,1,iF]) / 
#         (Global.Metric.dXdxI[iP,jP,1,1,3,3,iF] + Global.Metric.dXdxI[iP,jP,2,1,3,3,iF])
        K[iP,jP,1] = CE * uStar[iP,jP] * Global.Grid.dzeta[1]
      end
    end  
    @inbounds for iz = nz : -1 : 1
      @inbounds for jP = 1 : OP
        @inbounds for iP = 1 : OP
          if p[iP,jP,iz] > p_pbl
            K[iP,jP,iz] = K[iP,jP,1]
          else
            K[iP,jP,iz] = K[iP,jP,1] * exp(-((p_pbl - p[iP,jP,iz]) / p_strato)^2)  
          end
        end
      end
    end
  else
    @. K = 1.0  
  end   
  @. K = K * Rho

end




