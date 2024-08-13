using KernelAbstractions

function JacSchurGPU!(J,U,CG,Metric,Phys,Cache,Global,Param,Equation::Models.EquationType)

  backend = get_backend(U)
  FT = eltype(U)

  Nz = size(U,1)
  NumG = size(U,2)

  NG = div(256,Nz)
  group = (Nz, NG)
  ndrange = (Nz, NumG)
  dPresdRhoTh = Global.Model.dPresdRhoTh
  @views p = Cache.AuxG[:,:,1]

  KJacSchurKernel! = JacSchurKernel!(backend,group)

  KJacSchurKernel!(dPresdRhoTh,J.JRhoW,J.JWRho,J.JWRhoTh,J.JRhoThW,U,p,Metric.dz,Phys,Param,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end


@kernel inbounds = true function JacSchurKernel!(dPresdRhoTh,JRhoW,JWRho,JWRhoTh,JRhoThW,@Const(U),@Const(p),@Const(dz),Phys,Param)
  iz, iC   = @index(Local, NTuple)
  Iz,IC = @index(Global, NTuple)

  NG = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]    

  if Iz < Nz && IC <= NumG
    RhoPos = 1
    ThPos = 5
    RhoL = U[Iz,IC,RhoPos]
    RhoR = U[Iz+1,IC,RhoPos]
    RhoThL = U[Iz,IC,ThPos] + p[Iz,IC]
    RhoThR = U[Iz+1,IC,ThPos] + p[Iz+1,IC]
    dzL = dz[Iz,IC]
    dzR = dz[Iz+1,IC]

    RhoF = (RhoL * dzL + RhoR * dzR) / (dzL + dzR)
    # JRhoW low bidiagonal matrix
    # First row diagonal
    # Second row lower diagonal
    JRhoW[1,Iz,IC] = -RhoF / dzL
    JRhoW[2,Iz,IC] = RhoF / dzR

    dPdThL = dPresdRhoTh(RhoThL)
    dPdThR = dPresdRhoTh(RhoThR)
    # JWRhoTh upper bidiagonal matrix
    # First row upper diagonal
    # Second row diagonal
    JWRhoTh[1,Iz,IC] = -dPdThR / RhoF / ( eltype(dz)(0.5) * (dzL + dzR))
    JWRhoTh[2,Iz,IC] = dPdThL / RhoF / ( eltype(dz)(0.5) * (dzL + dzR))

    # JWRho upper bidiagonal matrix
    # First row upper diagonal
    # Second row diagonal
    JWRho[1,Iz,IC] = -Phys.Grav * dzR / RhoF / (dzL + dzR)
    JWRho[2,Iz,IC] = -Phys.Grav * dzL / RhoF / (dzL + dzR)

    RhoThF = (RhoThL * dzL + RhoThR * dzR) / (dzL + dzR)
    # JRhoThW low bidiagonal matrix
    # First row diagonal
    # Second row lower diagonal
    JRhoThW[1,Iz,IC] = -RhoThF / dzL
    JRhoThW[2,Iz,IC] = RhoThF / dzR
  end
end

#@inline function dPresdThGPU(RhoTh, Phys)
#  dpdTh = Phys.Rd * (Phys.Rd * RhoTh / Phys.p0)^(Phys.kappa / (eltype(RhoTh)(1) - Phys.kappa))
#end

function JacSchur!(J,U,CG,Phys,Global,Param,::Val{:Conservative})
  (  RhoPos,
      uPos,
      vPos,
      wPos,
      ThPos,
      NumV,
      NumTr) = Global.Model
  nz=Global.Grid.nz
  NF=Global.Grid.NumFaces
  nCol=size(U,2)
  nJ=nCol*nz

  D = J.CacheCol2
  Dp = J.CacheCol2
  Dm = J.CacheCol3
  dPdTh = J.CacheCol1
  K = J.CacheCol1

  for iC=1:nCol
    @views Pres = Global.Cache.AuxG[:,iC,1]
    @views Rho = U[:,iC,RhoPos]
    @views Th = U[:,iC,ThPos]
    @views Tr = U[:,iC,NumV+1:end]
    @views dz = Global.Metric.dz[:,iC]

    @views @. J.JRhoW[1,:,iC] = -1.0 / dz
    @views @. J.JRhoW[2,1:nz-1,iC] = 1.0 / dz[2:nz]

    dPresdTh!(dPdTh, Th, Rho, Tr, Phys, Global)
    @views @. Dp[1:nz-1] = dPdTh[2:nz] /  (1/2 * (dz[1:nz-1] + dz[2:nz]))
    @views @. Dm[1:nz-1] = dPdTh[1:nz-1] / (1/2 * (dz[1:nz-1] + dz[2:nz]))
    @views @. J.JWTh[1,2:nz,iC] = -Dp[1:nz-1]
    @views @. J.JWTh[2,:,iC] = Dm

    if Global.Model.Equation == "CompressibleMoist"
      dPresdRhoV!(dPdTh, Th, Rho, Tr, Pres, Phys, Global)
      @views @. Dp[1:nz-1] = dPdTh[2:nz] /  (1/2 * (dz[1:nz-1] + dz[2:nz]))
      @views @. Dm[1:nz-1] = dPdTh[1:nz-1] / (1/2 *(dz[1:nz-1] + dz[2:nz]))
      @views @. J.JWRhoV[1,2:nz,iC] = -Dp[1:nz-1]
      @views @. J.JWRhoV[2,:,iC] = Dm
    end  

    @views @. D[1:nz-1] = 1/2 * Phys.Grav 
    @views @. J.JWRho[1,2:nz,iC] = -D[1:nz-1] 
    @views @. J.JWRho[2,1:nz,iC] = -D 

    @views @. D[1:nz-1] = -(Th[1:nz-1] / Rho[1:nz-1] * dz[1:nz-1] + Th[2:nz] / Rho[2:nz] * dz[2:nz]) /
      (dz[1:nz-1] + dz[2:nz])
    @views @. J.JThW[1,:,iC] = D / dz
    @views @. J.JThW[2,1:nz-1,iC] = -D[1:nz-1] / dz[2:nz]
    if Global.Model.Thermo == "TotalEnergy" 
      @views @. D[1:nz-1] -= 1/2*(Pres[1:nz-1] + Pres[2:nz]) 
    elseif Global.Model.Thermo == "InternalEnergy"
      @views @. D = Pres / dz
      @views @. J.JThW[1,:,iC] = J.JThW[1,:,iC] - D
      @views @. J.JThW[2,1:nz-1,iC] = J.JThW[2,1:nz-1,iC] + D[2:nz]
    end  

    for iT = 1 : NumTr
      @views @. D[1:nz-1] = -(Tr[1:nz-1,iT] / Rho[1:nz-1] * dz[1:nz-1] + Tr[2:nz,iT] / Rho[2:nz] * dz[2:nz]) /
        (dz[1:nz-1] + dz[2:nz])
      @views @. J.JTrW[1,:,iC,iT] = D / dz
      @views @. J.JTrW[2,1:nz-1,iC,iT] = -D[1:nz-1] / dz[2:nz]
    end    

    if Global.Model.Damping
      @views DampingKoeff!(J.JWW[1,:,iC],CG,Global)
    end
    if Global.Model.VerticalDiffusion
      @views JDiff = J.JDiff[:,:,iC]
      @views KV = Global.Cache.AuxG[:,iC,2]
      # The Rho factor is already included in KV
      @views @. DF = - (KV[1:nz-1] + KV[2:nz]) / (dz[1:nz-1] + dz[2:nz])
      # J tridiagonal matrix
      # First row upper diagonal
      # Second row diagonal
      # Third row lower diagonal
      @views @. JDiff[1,1:nz-1] = DF / Rho[2:nz] / dz[1:nz-1]
      @views @. JDiff[3,2:nz] = DF / Rho[1:nz-1] / dz[2:nz]
      @views @. JDiff[2,:] = JDiff[1,:] + JDiff[3,:]
    end    
  end
end
