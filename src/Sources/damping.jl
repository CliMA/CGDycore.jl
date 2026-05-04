abstract type DampingValue end

Base.@kwdef struct DampingW <: DampingValue end

function (::DampingW)(H,StrideDamp,Relax,uPos,vPos,wPos,::Examples.VelocityS,::Grids.SphericalGrid)
  @inline function Damping(F,U,xS,z)
    FT = eltype(z)
    if z>=H-StrideDamp
      Damp = Relax *
        sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
      F[wPos] += -Damp * U[wPos]
    end
  end
  return Damping
end

function (::DampingW)(H,StrideDamp,Relax,uPos,vPos,wPos,::Examples.VelocityC,::Grids.SphericalGrid)
  @inline function Damping(F,U,xS,z)
    FT = eltype(F)
    if z>=H-StrideDamp
      Damp = Relax *
        sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
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
      wS = rot31 * uC + rot32 * vC + rot33 * wC
      FwS = -Damp * wS 
      FuC = rot31 * FwS
      FvC = rot32 * FwS
      FwC = rot33 * FwS
      F[uPos] += FuC
      F[vPos] += FvC
      F[wPos] += FwC
    end
  end
  return Damping
end

function Damping!(Damp,F,U,Aux,FE,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(F,1)
  Nz = size(F,2)
  NumG = size(F,3)
  NumV = size(F,4)
  NumGG = min(div(NumberThreadGPU,Nz*M),NumG)
  group = (M,Nz,NumGG)
  ndrange = (M,Nz,NumG)
  KDampKernel! = DampKernel!(backend, group)
  KDampKernel!(Damp,F,U,Aux,Metric.xS,Metric.zP,Val(NumV);ndrange=ndrange)
end

@kernel inbounds = true function DampKernel!(Damp,F,@Const(U),@Const(Aux),@Const(xS),@Const(zP),
  ::Val{NUMV}) where {NUMV}

  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)
  ULoc = @private eltype(F) (NUMV,)
  xSLoc = @private eltype(F) (2,)

  if ID <= ND
    zPLoc = zP[Iz,ID]
    xSLoc[1] = xS[1,ID]
    xSLoc[2] = xS[2,ID]
    @unroll for iv = 1 : NUMV
      ULoc[iv] = U[K,Iz,ID,iv]
      FLoc[iv] = eltype(F)(0)
    end
    Damp(FLoc,ULoc,xSLoc,zPLoc)
    @unroll for iv = 1 : NUMV
      F[K,Iz,ID,iv] += FLoc[iv]
    end
  end
end

@kernel inbounds = true function DampCGKernel!(Damp,F,@Const(U),@Const(zP))
  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    @views Fu,Fv,Fw = Damp(zP[Iz,ID],U[K,Iz,ID,:])
    F[K,Iz,ID,2] += Fu
    F[K,Iz,ID,3] += Fv
    F[K,Iz,ID,4] += Fw
  end
end

function Damping!(Damp,F,U,FE::FiniteElements.DGElement,Metric,NumberThreadGPU)
  backend = get_backend(F)
  M = size(F,1)
  Nz = size(F,2)
  DoF = size(FE.Glob,1)
  NF = size(FE.Glob,2)
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KDampKernel! = DampDGKernel!(backend, group)
  KDampKernel!(Damp,F,U,Metric.X,FE.Glob;ndrange=ndrange)
end

@kernel inbounds = true function DampDGKernel!(Damp,F,@Const(U),@Const(X),@Const(Glob))
  _,_,iD,  = @index(Local, NTuple)
  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    ind = Glob[ID,IF]
    @views Fu,Fv,Fw = Damp(X[ID,K,:,Iz,IF],U[K,Iz,ind,:])
    F[K,Iz,ind,2] += Fu
    F[K,Iz,ind,3] += Fv
    F[K,Iz,ind,4] += Fw
  end
end


