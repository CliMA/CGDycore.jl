abstract type DampingValue end

Base.@kwdef struct DampingW <: DampingValue end


# Variable positions as Int or Val
@inline PosVal(::Val{P}) where {P} = P
@inline PosVal(P::Integer) = P

# RelaxH [1/s]: optional sponge for the horizontal wind with the same vertical profile
function (::DampingW)(H,StrideDamp,Relax,uPosV,vPosV,wPosV,
  ::Union{Examples.VelocityC,Examples.VelocityS},::Grids.CartesianGrid;RelaxH=zero(Relax))
  uPos = PosVal(uPosV)
  vPos = PosVal(vPosV)
  wPos = PosVal(wPosV)
  @inline function Damping(F,U,xS,z)
    FT = eltype(z)
    if z>=H-StrideDamp
      Prof = sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
      F[wPos] += -Relax * Prof * U[wPos]
      if RelaxH > 0
        F[uPos] += -RelaxH * Prof * U[uPos]
        F[vPos] += -RelaxH * Prof * U[vPos]
      end
    end
  end
  return Damping
end

function (::DampingW)(H,StrideDamp,Relax,uPosV,vPosV,wPosV,::Examples.VelocityS,::Grids.SphericalGrid;
  RelaxH=zero(Relax))
  uPos = PosVal(uPosV)
  vPos = PosVal(vPosV)
  wPos = PosVal(wPosV)
  @inline function Damping(F,U,xS,z)
    FT = eltype(z)
    if z>=H-StrideDamp
      Prof = sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
      F[wPos] += -Relax * Prof * U[wPos]
      if RelaxH > 0
        F[uPos] += -RelaxH * Prof * U[uPos]
        F[vPos] += -RelaxH * Prof * U[vPos]
      end
    end
  end
  return Damping
end

function (::DampingW)(H,StrideDamp,Relax,uPosV,vPosV,wPosV,::Examples.VelocityC,::Grids.SphericalGrid;
  RelaxH=zero(Relax))
  uPos = PosVal(uPosV)
  vPos = PosVal(vPosV)
  wPos = PosVal(wPosV)
  @inline function Damping(F,U,xS,z)
    FT = eltype(F)
    if z>=H-StrideDamp
      Prof = sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
      Damp = Relax * Prof
      DampH = RelaxH * Prof
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
      # radial part with Relax, horizontal part (U - wS*e_r) with RelaxH
      FwS = -Damp * wS + DampH * wS
      FuC = rot31 * FwS - DampH * uC
      FvC = rot32 * FwS - DampH * vC
      FwC = rot33 * FwS - DampH * wC
      F[uPos] += FuC
      F[vPos] += FvC
      F[wPos] += FwC
    end
  end
  return Damping
end

# IzStart: lowest level that reaches the damping layer in any column, the levels below are skipped
function Damping!(Damp,F,U,Aux,FE,Metric,NumberThreadGPU;IzStart=1)
  backend = get_backend(F)
  M = size(F,1)
  Nz = size(F,2) - IzStart + 1
  NumG = size(F,3)
  NumV = size(F,4)
  NumGG = min(div(NumberThreadGPU,Nz*M),NumG)
  group = (M,Nz,NumGG)
  ndrange = (M,Nz,NumG)
  KDampKernel! = DampKernel!(backend, group)
  KDampKernel!(Damp,F,U,Aux,Metric.xS,Metric.zP,Val(NumV),IzStart - 1;ndrange=ndrange)
end

# lowest level with zP >= H - StrideDamp in any column
function DampingIzStart(zP,H,StrideDamp)
  zPH = Array(zP)
  for Iz in axes(zPH,1)
    if maximum(view(zPH,Iz,:)) >= H - StrideDamp
      return Iz
    end
  end
  return size(zPH,1)
end

@kernel inbounds = true function DampKernel!(Damp,F,@Const(U),@Const(Aux),@Const(xS),@Const(zP),
  ::Val{NUMV},IzOffset) where {NUMV}

  _,_,iD  = @index(Local, NTuple)
  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  FLoc = @private eltype(F) (NUMV,)
  xSLoc = @private eltype(F) (2,)

  Iz += IzOffset
  if ID <= ND
    zPLoc = zP[Iz,ID]
    xSLoc[1] = xS[1,ID]
    xSLoc[2] = xS[2,ID]
    @unroll for iv = 1 : NUMV
      FLoc[iv] = eltype(F)(0)
    end
    # U as view: only the components the damping uses are loaded, only changed ones written
    Damp(FLoc,view(U,K,Iz,ID,:),xSLoc,zPLoc)
    @unroll for iv = 1 : NUMV
      if FLoc[iv] != eltype(F)(0)
        F[K,Iz,ID,iv] += FLoc[iv]
      end
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


