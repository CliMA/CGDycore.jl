function held_suarez_physics(u::Real, v::Real, T::Real, p::Real, lat::Real)
    # Constants
    p0 = 100000.0
    day_to_sec = 86400.0
    sigma_b = 0.7

    # 1. Coordinate setup
    sigma = p / p0
    cos_lat = cos(lat)
    sin_lat = sin(lat)

    # 2. Temperature Forcing (Newtonian Relaxation)
    ka = 1.0 / (40.0 * day_to_sec)
    ks = 1.0 / (4.0 * day_to_sec)

    # Te profile calculation
    term1 = 315.0 - 60.0 * sin_lat^2
    term2 = 10.0 * log(sigma) * cos_lat^2
    Te = max(200.0, (term1 - term2) * sigma^(2/7))

    # Temperature relaxation rate
    kT = ka + (ks - ka) * cos_lat^4 * max(0.0, (sigma - sigma_b) / (1.0 - sigma_b))
    dT_dt = -kT * (T - Te)

    # 3. Wind Forcing (Rayleigh Friction)
    # Friction only acts near the surface (sigma > sigma_b)
    kf = 1.0 / (1.0 * day_to_sec) # 1 day timescale at surface
    kv = kf * max(0.0, (sigma - sigma_b) / (1.0 - sigma_b))

    du_dt = -kv * u
    dv_dt = -kv * v

    return du_dt, dv_dt, dT_dt
end

abstract type ForcingType end

Base.@kwdef struct HeldSuarezDryForcing <: ForcingType end

function (::HeldSuarezDryForcing)(Param,RhoPos,uPos,vPos,wPos,ThPos,pPos,::Examples.VelocityS)
  @inline function Force(F,U,Aux,xS)
    FT = eltype(U)
    p = Aux[pPos]
    lat = xS[2]
    Sigma = p / FT(P.p0)
    SigmaPowKappa = fast_powGPU(Sigma,FT(P.kappa))
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    sinlat = sin(lat)
    coslat = cos(lat)
    F[uPos] += -(Param.k_f * height_factor) * U[uPos]
    F[vPos] += -(Param.k_f * height_factor) * U[vPos]
    kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * coslat * coslat * coslat * coslat
    Teq = (Param.T_equator - Param.DeltaT_y * sinlat * sinlat -
      Param.DeltaTh_z * log(Sigma) * coslat * coslat) * SigmaPowKappa
    Teq = max(Param.T_min, Teq)
    T = p / (U[RhoPos] * FT(P.Rd))
    DeltaT =  kT * (T - Teq)
    F[ThPos]  += - U[ThPos] / T * DeltaT

#=
    u = U[uPos]
    v = U[vPos]
    p0 = 100000.0
    day_to_sec = 86400.0
    sigma_b = 0.7

    # 1. Coordinate setup
    sigma = p / p0
    cos_lat = cos(lat)
    sin_lat = sin(lat)

    # 2. Temperature Forcing (Newtonian Relaxation)
    ka = 1.0 / (40.0 * day_to_sec)
    ks = 1.0 / (4.0 * day_to_sec)

    # Te profile calculation
    term1 = 315.0 - 60.0 * sin_lat^2
    term2 = 10.0 * log(sigma) * cos_lat^2
    Te = max(200.0, (term1 - term2) * sigma^(2/7))

    # Temperature relaxation rate
    kT = ka + (ks - ka) * cos_lat^4 * max(0.0, (sigma - sigma_b) / (1.0 - sigma_b))
    dT_dt = -kT * (T - Te)

    # 3. Wind Forcing (Rayleigh Friction)
    # Friction only acts near the surface (sigma > sigma_b)
    kf = 1.0 / (1.0 * day_to_sec) # 1 day timescale at surface
    kv = kf * max(0.0, (sigma - sigma_b) / (1.0 - sigma_b))

    du_dt = -kv * u
    dv_dt = -kv * v
    if abs(du_dt-F[uPos])>0.001
      @show du_dt,F[uPos]
      stop
    end  
    if abs(dv_dt-F[vPos])>0.001
      @show dv_dt,F[vPos]
      stop
    end  
    if abs(dT_dt-DeltaT)>0.001
      @show dv_dt,DeltaT
      stop
    end  
    =#
  end
  return Force
end

function (::HeldSuarezDryForcing)(Param,RhoPos,uPos,vPos,wPos,ThPos,pPos,::Examples.VelocityC)
  @inline function Force(F,U,Aux,xS)
    FT = eltype(U)
    p = Aux[pPos]
    lon = xS[1]
    lat = xS[2]
    sinlat = sin(lat)
    coslat = cos(lat)
    sinlon = sin(lon)
    coslon = cos(lon)
    rot11 = -sinlon
    rot21 = -sinlat * coslon
    rot31 = coslat * coslon
    rot12 = coslon
    rot22 = -sinlat * sinlon
    rot32 = coslat * sinlon
    rot13 = eltype(lon)(0)
    rot23 = coslat
    rot33 = sinlat
    uC = U[uPos]
    vC = U[vPos]
    wC = U[wPos]
    uS = rot11 * uC + rot12 * vC + rot13 * wC
    vS = rot21 * uC + rot22 * vC + rot23 * wC
    Sigma = p / FT(P.p0)
    SigmaPowKappa = fast_powGPU(Sigma,FT(P.kappa))
    height_factor = max(FT(0), (Sigma - Param.sigma_b) / (FT(1) - Param.sigma_b))
    FuS = -(Param.k_f * height_factor) * uS
    FvS = -(Param.k_f * height_factor) * vS
    FuC = rot11 * FuS + rot21 * FvS
    FvC = rot12 * FuS + rot22 * FvS
    FwC = rot13 * FuS + rot23 * FvS
    F[uPos] += FuC
    F[vPos] += FvC
    F[wPos] += FwC
    kT = Param.k_a + (Param.k_s - Param.k_a) * height_factor * coslat * coslat * coslat * coslat
    Teq = (Param.T_equator - Param.DeltaT_y * sinlat * sinlat -
      Param.DeltaTh_z * log(Sigma) * coslat * coslat) * SigmaPowKappa
    Teq = max(Param.T_min, Teq)
    T = p / (U[RhoPos] * FT(P.Rd))
    DeltaT =  kT * (T - Teq)
    F[ThPos]  += - U[ThPos] / T * DeltaT
  end
  return Force
end

function Forcing!(Force,F,U,Aux,FE,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(F,1)
  Nz = size(F,2)
  NumG = size(F,3)
  NumV = size(F,4)
  NumAux = size(Aux,4)
  NumGG = min(div(NumberThreadGPU,Nz*M),NumG)
  group = (M,Nz,NumGG)
  ndrange = (M,Nz,NumG)
  KForceKernel! = ForceKernel!(backend, group)
  KForceKernel!(Force,F,U,Aux,Metric.xS,Val(NumV),Val(NumAux);ndrange=ndrange)
end

@kernel inbounds = true function ForceKernel!(Force,F,@Const(U),@Const(Aux),@Const(xS),
  ::Val{NUMV}, ::Val{NUMAUX}) where {NUMV,NUMAUX}

  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)
  ULoc = @private eltype(F) (NUMV,)
  AuxLoc = @private eltype(F) (NUMAUX,)
  xSLoc = @private eltype(F) (2,)

  if ID <= ND
    xSLoc[1] = xS[1,ID]
    xSLoc[2] = xS[2,ID]
    @unroll for iv = 1 : NUMV
      ULoc[iv] = U[K,Iz,ID,iv]
      FLoc[iv] = eltype(F)(0)
    end  
    @unroll for iAux = 1 : NUMAUX
      AuxLoc[iAux] = Aux[K,Iz,ID,iAux]
    end  
    Force(FLoc,ULoc,AuxLoc,xSLoc)
    @unroll for iv = 1 : NUMV
      F[K,Iz,ID,iv] += FLoc[iv]
    end  
  end
end 

@inline fast_powGPU(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))
