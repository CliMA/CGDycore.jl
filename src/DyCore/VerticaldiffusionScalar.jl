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

function BoundaryFluxSaclar!(Fc,c,cS,CG,Global,iF)
  @views @. Fc -= Global.Model.Param.CTr * Global.Cache.uStar[:,:,iF] * (c - cS) *
    Global.Metric.dXdxIF[:,:,1,3,3,iF]
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



