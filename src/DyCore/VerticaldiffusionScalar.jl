function VerticalDiffusionSaclar!(Fc,c,Rho,K,CG,Global,iF)

nz = Global.Grid.nz
gradc = Global.Cache.CacheC1
@views dXdxIF33 = Global.Metric.dXdxIF[:,:,:,3,3,iF] 
@views JF = Global.Metric.JF[:,:,:,iF] 
qCG = Global.Cache.CacheC3
@. qCG = c / Rho 
# Gradient
  @inbounds for iz = 1:nz-1
    @views @. gradc[:,:,iz+1] = 0.5 * (K[:,:,iz] + K[:,:,iz+1]) * (c[:,:,iz+1] - c[:,:,iz]) * dXdxIF33[:,:,iz+1] / JF[:,:,iz+1] 
  end
# Divergence  
  @views @. Fc[:,:,1] += 0.5 * gradc[:,:,2] * dXdxIF33[:,:,2]
  @inbounds for iz = 2:nz-1
    @views @. Fc[:,:,iz] += 0.5 * (gradc[:,:,iz+1] * dXdxIF33[:,:,iz+1] - gradc[:,:,iz] * dXdxIF33[:,:,iz])
  end
  @views @. Fc[:,:,nz] -= 0.5 * gradc[:,:,nz] * dXdxIF33[:,:,nz]
end

function BoundaryFluxScalar!(Fc,c,cS,CG,Global,iF)
  @views @. Fc -= Global.Model.Param.CTr * Global.Cache.uStar[:,:,iF] * (c - cS) *
    Global.Metric.dXdxIF[:,:,1,3,3,iF]
end

function BoundaryFluxScalar!(Fc,Th,Rho,Tr,CG,Global,iF)
  if Global.Model == "HeldSuarezMoist"
    println(" Hello ")  
    ThPos=Global.Model.ThPos  
    RhoPos=Global.Model.RhoPos  
    RhoVPos=Global.Model.RhoVPos  
    NumV=Global.Model.NumV  
    TSurf = Global.Cache.TSurf
    uStar = Global.Cache.uStar
    CE = Global.Model.Param.CE
    CH = Global.Model.Param.CH
    Rd = Global.Phys.Rd
    Cpd = Global.Phys.Cpd
    Rv = Global.Phys.Rv
    Cpv = Global.Phys.Cpv
    p0 = Global.Phys.p0
    p = Global.Cache.Pres(:,:,1,iF)
    T = Global.Cache.Temp(:,:,1,iF)
    dXdxIF = Global.Metric.dXdxIF[:,:,1,iF,3,3]
    JF = Global.Metric.JF[:,:,1,iF]
    @. BoundaryFluxHeldSuarez!(Fc[:,:,ThPos],Fc[:,:,RhoPos],Fc[:,:,RhoVPos+NumV],
      Th,Rho,Tr[:,:,RhoVPos],p,T,dXdxIF,JF,CH,CE,TSurf ,uStar, Rd, Cpd, Rv, Cpv, p0)
    stop
      
  end    
end

BoundaryFluxHeldSuarez!(FTh,FRho,FRhoV,Th,Rho,RhoV,TSea,p,T,dXdxIF,JF,CH,CE,TS,uStar, Rd, Cpd, Rv, Cpv, p0)
  T_C = TS - 273.15
  p_vs = 611.2 * exp(17.62 * T_C / (243.12 + T_C))
  RhoVSurface = Rho * p_vs / (Rv * TS) 
  LatFlux = 0.25 * CE * uStar * dXdxIF^2 /JF  *(RhoV(ix,iy,iz,1)-RhoVSurface) 
  SensFlux = 0.25 * CH * uStar * dXdxIF^2 /JF  *(T-Ts) 
  FRho = LatFlux
  FRhoV = LatFlux
  RhoDry = Rho - RhoV
  rrv=RhoV/RhoDry
  Rm=Rd+rrv*Rv
  Cpml=Cpd+Cpv*rrv
  PrePi=(p / p0)^(Rm / Cpml)
  FTh = Th / T * SensFlux + ((Rv / Rm) - log(PrePi)*(Rv / Rm - Cpv / Cpml)) / Rho * MoistFlux  
end


function uStarCoefficient!(uStar,U,V,WC,CG,Global,iF)
# Computation norm_v_a
# |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)
  @views @. uStar =  sqrt(U * U + V * V + WC * WC -
   (Global.Metric.nS[:,:,1,iF] * U + Global.Metric.nS[:,:,2,iF] * V + Global.Metric.nS[:,:,3,iF] * WC)^2)
end

function eddy_diffusivity_coefficient!(K,U,V,WC,Rho,CG,Global,iF) 
  if Global.Model.Problem == "HeldSuarezMoist"
    C_E = Global.Model.Param.C_E 
    p_pbl = Global.Model.Param.p_pbl 
    p_strato = Global.Model.Param.p_strato 
#   Computation norm_v_a  
#   |v_a| = |v - n(n*v)| = sqrt(v*v -(n*v)^2)  
    @views uStar = Global.Cache.uStar[:,:,iF]
    @views @. K[:,:,1] = 0.5 * C_E * uStar * Global.Metric.dz[:,:,1,iF]
    @views p = Global.Cache.Pres[:,:,:,iF]
    for iz = size(K,3) : -1 : 1
      for jP = 1 : size(K,2)
        for iP = 1 : size(K,1)
          if p[iP,jP,iz] < p_pbl
            K[iP,jP,iz] = K[iP,jP,1]
          else
            K[iP,jP,iz] = K[iP,jP,1] * exp(-((p_pbl - p[iP,jP,iz]) / p_strato)^2)  
          end
        end
      end
    end
  else
    @views @. K = 1.0  
  end   
  @. K = K * Rho
end



