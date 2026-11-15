# Kernels for the Cartesian velocity frame (VelocityFrame Cartesian): the tangential
# velocity is stored in Cartesian components (ux,uy,uz) at (uPos,vPos,uzPos), the radial
# velocity w at wPos. Each row J grad xi^a of the Cartesian metric is split into its
# tangential part dXdxI[a,1:3] and its normal part dXdxI[a,4] = J grad xi^a . n, n = X/|X|,
# so that the contravariant velocity is
#   J grad xi^a . (uh + w n) = dXdxI[a,1:3] . uh + dXdxI[a,4] w,
# and gradients add only tangential vectors to uh (no cancellation in Float32).
# No rotation matrix and no trigonometric functions are needed.

# Options that still use the spherical components u, v (surface fluxes, momentum diffusion,
# forcing, turbulence, EDMF, limiter, energy states) are not available in the Cartesian frame
CheckVelocityFrame(Model,::Models.SphericalFrame) = nothing

function CheckVelocityFrame(Model,::Models.CartesianFrame)
  Model.State == "Dry" || throw(ArgumentError("VelocityFrame Cartesian: only State Dry"))
  for (flag, name) in ((Model.SurfaceFlux, "SurfaceFlux"), (Model.SurfaceFluxMom, "SurfaceFluxMom"),
    (Model.VerticalDiffusion, "VerticalDiffusion"), (Model.VerticalDiffusionMom, "VerticalDiffusionMom"),
    (Model.Forcing, "Forcing"), (Model.Turbulence, "Turbulence"), (Model.EDMF, "EDMF"),
    (Model.HorLimit, "HorLimit"))
    flag && throw(ArgumentError("VelocityFrame Cartesian: $name is not available"))
  end
  return nothing
end

# horizontal contravariant velocities, both vertical nodes of the layer (w drops out,
# J grad xi^a . n = 0 for a = 1,2 since the zeta lines are radial)
@inline function Contra12Cart(Rho,u1,u2,u3,dXdxI)
  uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u1 +
    (dXdxI[1,2,1] + dXdxI[1,2,2]) * u2 +
    (dXdxI[1,3,1] + dXdxI[1,3,2]) * u3)
  vCon = Rho * ((dXdxI[2,1,1] + dXdxI[2,1,2]) * u1 +
    (dXdxI[2,2,1] + dXdxI[2,2,2]) * u2 +
    (dXdxI[2,3,1] + dXdxI[2,3,2]) * u3)
  return uCon, vCon
end

# vertical contravariant mass flux at the interface between two layers,
# dXdxI = dXdxI[3,:,:,ID,Iz:Iz+1,IF]
@inline function Contra3Cart(Rho,u1,u2,u3,w,dXdxI)
  wCon = Rho[1] * (dXdxI[1,2,1] * u1[1] + dXdxI[2,2,1] * u2[1] + dXdxI[3,2,1] * u3[1] +
    dXdxI[4,2,1] * w) +
    Rho[2] * (dXdxI[1,1,2] * u1[2] + dXdxI[2,1,2] * u2[2] + dXdxI[3,1,2] * u3[2] +
    dXdxI[4,1,2] * w)
end

# DivRhoTrUpwind3New2Kernel! in the Cartesian frame
@kernel inbounds = true function DivRhoTrUpwind3New2CartKernel!(F,NumV,NumTr,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),uzPos)

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  uzCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  JCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  MCCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  if Iz <= Nz
    cCol[I,J,iz+1] = U[1,Iz,ind,5] / U[1,Iz,ind,1]
    JCol[I,J,iz+1] = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    RhoCol[I,J,iz] = U[1,Iz,ind,1]
    uCol[I,J,iz] = U[1,Iz,ind,2]
    vCol[I,J,iz] = U[1,Iz,ind,3]
    uzCol[I,J,iz] = U[1,Iz,ind,uzPos]
    wCol[I,J,iz] = U[1,Iz,ind,4]
    MCCol[I,J,iz] = M[Iz,ind,1] + M[Iz,ind,2]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = U[1,Izm1,ind,5] / U[1,Izm1,ind,1]
    JCol[I,J,iz] = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = U[1,Izp1,ind,5] / U[1,Izp1,ind,1]
    JCol[I,J,iz+2] = JJ[ID,1,Izp1,IF] + JJ[ID,2,Izp1,IF]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = U[1,Izp2,ind,5] / U[1,Izp2,ind,1]
    JCol[I,J,iz+3] = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    RhoCol[I,J,iz+1] = U[1,Izp1,ind,1]
    uCol[I,J,iz+1] = U[1,Izp1,ind,2]
    vCol[I,J,iz+1] = U[1,Izp1,ind,3]
    uzCol[I,J,iz+1] = U[1,Izp1,ind,uzPos]
    MCCol[I,J,iz+1] = M[Izp1,ind,1] + M[Izp1,ind,2]
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz < Nz
    @views wCon = Contra3Cart(RhoCol[I,J,iz:iz+1],uCol[I,J,iz:iz+1],vCol[I,J,iz:iz+1],
      uzCol[I,J,iz:iz+1],wCol[I,J,iz],dXdxI[3,:,:,ID,Iz:Iz+1,IF])
    wCol[I,J,iz] = wCon
    cFL, cFR = RecU4(cCol[I,J,iz],cCol[I,J,iz+1],cCol[I,J,iz+2],cCol[I,J,iz+3],
      JCol[I,J,iz],JCol[I,J,iz+1],JCol[I,J,iz+2],JCol[I,J,iz+3])
#   central second order at the interfaces next to the boundaries, no ghost cells
    if Iz == 1 || Iz == Nz - 1
      cFL = eltype(F)(0.5) * (cCol[I,J,iz+1] + cCol[I,J,iz+2])
      cFR = cFL
    end
    Flux = eltype(F)(0.25) * ((abs(wCon) + wCon) * cFL +
      (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic F[1,Iz,ind,5] += -Flux / MCCol[I,J,iz]
    @atomic :monotonic F[1,Iz+1,ind,5] += Flux / MCCol[I,J,iz+1]
    Flux = eltype(F)(0.5)*wCon
    @atomic :monotonic F[1,Iz,ind,1] += -Flux / MCCol[I,J,iz]
    @atomic :monotonic F[1,Iz+1,ind,1] += Flux / MCCol[I,J,iz+1]
  end

  if Iz <= Nz
    uCol[I,J,iz], vCol[I,J,iz] = Contra12Cart(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],uzCol[I,J,iz],
      view(dXdxI,1:2,1:3,:,ID,Iz,IF))
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]
  if Iz <= Nz
    DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
    DivTr1 = D[I,1] * cCol[1,J,iz+1]
    DivTr2 = D[J,1] * cCol[I,1,iz+1]
    DivRho = D[I,1] * uCol[1,J,iz] + D[J,1] * vCol[I,1,iz]
    for k = 2 : N
      DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
      DivTr1 += D[I,k] * cCol[k,J,iz+1]
      DivTr2 += D[J,k] * cCol[I,k,iz+1]
      DivRho += D[I,k] * uCol[k,J,iz] + D[J,k] * vCol[I,k,iz]
    end
    @atomic :monotonic F[1,Iz,ind,5] += eltype(F)(0.5) * (DivRhoTr + uCol[I,J,iz] * DivTr1 + vCol[I,J,iz] * DivTr2 +
      cCol[I,J,iz+1] * DivRho)/ MCCol[I,J,iz]
    @atomic :monotonic F[1,Iz,ind,1] += DivRho / MCCol[I,J,iz]
  end

  for iT = NumV + 1 : NumV + NumTr
#   Second tracer
    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    if Iz <= Nz
      cCol[I,J,iz+1] = U[1,Iz,ind,iT] / U[1,Iz,ind,1]
    end
    if iz == 1
      Izm1 = max(Iz - 1,1)
      cCol[I,J,iz] = U[1,Izm1,ind,iT] / U[1,Izm1,ind,1]
    end
    if iz == ColumnTilesDim || Iz == Nz
      Izp1 = min(Iz + 1,Nz)
      cCol[I,J,iz+2] = U[1,Izp1,ind,iT] / U[1,Izp1,ind,1]
      Izp2 = min(Iz + 2,Nz)
      cCol[I,J,iz+3] = U[1,Izp2,ind,iT] / U[1,Izp2,ind,1]
    end
    @synchronize

    ID = I + (J - 1) * N
    ind = Glob[ID,IF]
    if Iz < Nz
      wCon = wCol[I,J,iz]
      cFL, cFR = RecU4(cCol[I,J,iz],cCol[I,J,iz+1],cCol[I,J,iz+2],cCol[I,J,iz+3],
        JCol[I,J,iz],JCol[I,J,iz+1],JCol[I,J,iz+2],JCol[I,J,iz+3])
      Flux = eltype(F)(0.25) * ((abs(wCon) + wCon) * cFL +
        (-abs(wCon) + wCon) * cFR)
      @atomic :monotonic F[1,Iz,ind,iT] += -Flux / MCCol[I,J,iz]
      @atomic :monotonic F[1,Iz+1,ind,iT] += Flux / MCCol[I,J,iz+1]
    end
    if Iz <= Nz
      DivRhoTr = D[I,1] * uCol[1,J,iz] * cCol[1,J,iz+1] + D[J,1] * vCol[I,1,iz] * cCol[I,1,iz+1]
      for k = 2 : N
        DivRhoTr += D[I,k] * uCol[k,J,iz] * cCol[k,J,iz+1] + D[J,k] * vCol[I,k,iz] * cCol[I,k,iz+1]
      end
      @atomic :monotonic F[1,Iz,ind,iT] += DivRhoTr / MCCol[I,J,iz]
    end
  end
end

# tangential part of row a of the Cartesian metric at vertical node K
@inline MetricRow(dXdxI,a,K,ID,Iz,IF) =
  SVector(dXdxI[a,1,K,ID,Iz,IF],dXdxI[a,2,K,ID,Iz,IF],dXdxI[a,3,K,ID,Iz,IF])
# full row J grad xi^a = tangential part + normal part n, for the cross products of the curl
@inline MetricRowFull(dXdxI,a,K,ID,Iz,IF,n) = MetricRow(dXdxI,a,K,ID,Iz,IF) + dXdxI[a,4,K,ID,Iz,IF] * n

# Lower boundary: no flow through the orography, J grad zeta . (uh + w n) = 0
@inline BottomW(dXdxI,u1,u2,u3,ID,IF) =
  -(dXdxI[3,1,1,ID,1,IF] * u1 + dXdxI[3,2,1,ID,1,IF] * u2 + dXdxI[3,3,1,ID,1,IF] * u3) /
  dXdxI[3,4,1,ID,1,IF]

# GradFullKernel! (Split = false) and GradFullSplitKernel! (Split = true) in the Cartesian
# frame. The tangential tendency is added as a Cartesian vector and projected after the
# assembly (ProjectTangential!), the radial component enters w with dXdxI[3,4].
@kernel inbounds = true function GradFullCartKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),@Const(X),
  @Const(JJ),@Const(M),@Const(Glob),Gravitation,uzPos,::Val{Split},Rd,Cpd) where {Split}

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  KinF = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+2)
  pCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  TCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  ThCol = @localmem eltype(F) (N,N,ColumnTilesDim)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    pCol[I,J,iz] = p[1,Iz,ind]  
    RhoCol[I,J,iz+1] = U[1,Iz,ind,1]
    if Split
      TCol[I,J,iz] = p[1,Iz,ind] / (Rd * U[1,Iz,ind,1])
      ThCol[I,J,iz] = U[1,Iz,ind,5] / U[1,Iz,ind,1]
    end
    KinF[I,J,1,iz+1] = eltype(F)(0.5) * (U[1,Iz,ind,2]^2 + U[1,Iz,ind,3]^2 + U[1,Iz,ind,uzPos]^2)
    KinF[I,J,2,iz+1] = KinF[I,J,1,iz+1] + eltype(F)(0.5) * U[1,Iz,ind,4]^2
    if iz == 1 && Iz == 1
      wCol = BottomW(dXdxI,U[1,Iz,ind,2],U[1,Iz,ind,3],U[1,Iz,ind,uzPos],ID,IF)
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * wCol^2
      KinF[I,J,2,iz] = KinF[I,J,1,iz+1]
    elseif iz == 1
      RhoCol[I,J,1] = U[1,Iz-1,ind,1]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[1,Iz-1,ind,4]^2
      KinF[I,J,2,iz] =  eltype(F)(0.5) * (U[1,Iz-1,ind,4]^2 + U[1,Iz-1,ind,2]^2 + U[1,Iz-1,ind,3]^2 +
        U[1,Iz-1,ind,uzPos]^2)
    else
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[1,Iz-1,ind,4]^2  
    end
    if iz == ColumnTilesDim && Iz < Nz
      pCol[I,J,iz+1] = p[1,Iz+1,ind]  
      RhoCol[I,J,iz+2] = U[1,Iz+1,ind,1]
      KinF[I,J,1,iz+2] = eltype(F)(0.5) * (U[1,Iz+1,ind,2]^2 + U[1,Iz+1,ind,3]^2 + U[1,Iz+1,ind,uzPos]^2)
      KinF[I,J,1,iz+2] +=  eltype(F)(0.5) * U[1,Iz,ind,4]^2
    end  
  end

  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    DXKinF1 = D[I,1] * KinF[1,J,1,iz+1]
    DYKinF1 = D[J,1] * KinF[I,1,1,iz+1]
    DXKinF2 = D[I,1] * KinF[1,J,2,iz+1]
    DYKinF2 = D[J,1] * KinF[I,1,2,iz+1]
    for k = 2 : N
      DXKinF1 += D[I,k] * KinF[k,J,1,iz+1]
      DYKinF1 += D[J,k] * KinF[I,k,1,iz+1]
      DXKinF2 += D[I,k] * KinF[k,J,2,iz+1]
      DYKinF2 += D[J,k] * KinF[I,k,2,iz+1]
    end
    Rho = RhoCol[I,J,iz+1]
#   horizontal pressure gradient in reference coordinates, divided by rho
    if Split
      DXT = D[I,1] * TCol[1,J,iz]
      DYT = D[J,1] * TCol[I,1,iz]
      DXPi = D[I,1] * TCol[1,J,iz] / ThCol[1,J,iz]
      DYPi = D[J,1] * TCol[I,1,iz] / ThCol[I,1,iz]
      DXTh = D[I,1] * ThCol[1,J,iz]
      DYTh = D[J,1] * ThCol[I,1,iz]
      for k = 2 : N
        DXT += D[I,k] * TCol[k,J,iz]
        DYT += D[J,k] * TCol[I,k,iz]
        DXPi += D[I,k] * TCol[k,J,iz] / ThCol[k,J,iz]
        DYPi += D[J,k] * TCol[I,k,iz] / ThCol[I,k,iz]
        DXTh += D[I,k] * ThCol[k,J,iz]
        DYTh += D[J,k] * ThCol[I,k,iz]
      end
      Th = ThCol[I,J,iz]
      Pi = TCol[I,J,iz] / Th
      GX = eltype(F)(0.5) * Cpd * (Th * DXPi + DXT - Pi * DXTh)
      GY = eltype(F)(0.5) * Cpd * (Th * DYPi + DYT - Pi * DYTh)
    else
      DXpCol = D[I,1] * pCol[1,J,iz]
      DYpCol = D[J,1] * pCol[I,1,iz]
      for k = 2 : N
        DXpCol += D[I,k] * pCol[k,J,iz]
        DYpCol += D[J,k] * pCol[I,k,iz]
      end
      GX = DXpCol / Rho
      GY = DYpCol / Rho
    end

    J11 = MetricRow(dXdxI,1,1,ID,Iz,IF)
    J12 = MetricRow(dXdxI,1,2,ID,Iz,IF)
    J21 = MetricRow(dXdxI,2,1,ID,Iz,IF)
    J22 = MetricRow(dXdxI,2,2,ID,Iz,IF)
    J31 = MetricRow(dXdxI,3,1,ID,Iz,IF)
    J32 = MetricRow(dXdxI,3,2,ID,Iz,IF)
    GradF1 = -(J11 * DXKinF1 + J21 * DYKinF1)
    GradF2 = -(J12 * DXKinF2 + J22 * DYKinF2)
    GradP = (J11 + J12) * GX + (J21 + J22) * GY

    GradwF1 = eltype(F)(0)
    GradwF2 = eltype(F)(0)
    if Iz > 1 
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+1] - KinF[I,J,2,iz] )  
      GradF1 += -GradZ * J31
      GradwF1 += -Rho * GradZ * dXdxI[3,4,1,ID,Iz,IF]
    end  
    if Iz < Nz
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+2] - KinF[I,J,2,iz+1] )  
      GradF2 += -GradZ * J32
      GradwF2 += -Rho * GradZ * dXdxI[3,4,2,ID,Iz,IF]
    end  
    GradZ = eltype(F)(0.5) * (KinF[I,J,2,iz+1] - KinF[I,J,1,iz+1])
    GradF2 += -GradZ * J32
    GradF1 += -GradZ * J31
    GradwF2 += -Rho * GradZ * dXdxI[3,4,2,ID,Iz,IF]
    GradwF1 += -Rho * GradZ * dXdxI[3,4,1,ID,Iz,IF]

    x = eltype(F)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(F)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(F)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    Grav = Gravitation(x,y,z)
    GradZ = -Grav * Rho *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,4,1,ID,Iz,IF] + dXdxI[3,4,2,ID,Iz,IF])
    GradG = GradZ * (J31 + J32)

    FH = (GradF1 + GradF2 - GradP - GradG / Rho) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,2] += FH[1]
    @atomic :monotonic F[1,Iz,ind,3] += FH[2]
    @atomic :monotonic F[1,Iz,ind,uzPos] += FH[3]

    if Iz < Nz
      GradZ = eltype(F)(0.5) * (pCol[I,J,iz+1] - pCol[I,J,iz])
      Gradw =  GradZ* (dXdxI[3,4,2,ID,Iz,IF] + dXdxI[3,4,1,ID,Iz+1,IF])
      x = eltype(F)(0.5) * (X[ID,2,1,Iz,IF] + X[ID,2,1,Iz+1,IF])
      y = eltype(F)(0.5) * (X[ID,2,2,Iz,IF] + X[ID,1,2,Iz+1,IF])
      z = eltype(F)(0.5) * (X[ID,2,3,Iz,IF] + X[ID,1,3,Iz+1,IF])
      Grav = Gravitation(x,y,z)
      Gradp = -(Gradw +
        Grav * (RhoCol[I,J,iz+1] * JJ[ID,2,Iz,IF] + RhoCol[I,J,iz+2] * JJ[ID,1,Iz+1,IF])) 
      @atomic :monotonic F[1,Iz,ind,4] += (GradwF2 + Gradp) / 
        (RhoCol[I,J,iz+2] * M[Iz+1,ind,1] + RhoCol[I,J,iz+1] * M[Iz,ind,2])
    end  
    if Iz > 1  
      @atomic :monotonic F[1,Iz-1,ind,4] += GradwF1 / 
        (RhoCol[I,J,iz] * M[Iz-1,ind,2] + RhoCol[I,J,iz+1] * M[Iz,ind,1])
    end  
  end

end

# Projection of the assembled tangential tendency onto the tangent plane,
# F_uh -> F_uh - (n . F_uh) n with the radial normal nR of the node
@kernel inbounds = true function ProjectTangentialKernel!(F,@Const(nR),uzPos)

  Iz,ind = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if Iz <= Nz && ind <= NumG
    n1 = nR[1,ind]
    n2 = nR[2,ind]
    n3 = nR[3,ind]
    s = n1 * F[1,Iz,ind,2] + n2 * F[1,Iz,ind,3] + n3 * F[1,Iz,ind,uzPos]
    F[1,Iz,ind,2] -= s * n1
    F[1,Iz,ind,3] -= s * n2
    F[1,Iz,ind,uzPos] -= s * n3
  end
end

function ProjectTangential!(F,Metric,Global,::Models.SphericalFrame)
end

function ProjectTangential!(F,Metric,Global,Frame::Models.CartesianFrame)
  backend = get_backend(F)
  Nz = size(F,2)
  NumG = size(F,3)
  NDoFG = min(div(Global.ParallelCom.NumberThreadGPU,Nz),NumG)
  KProjectTangentialKernel! = ProjectTangentialKernel!(backend,(Nz,NDoFG))
  KProjectTangentialKernel!(F,Metric.nR,Frame.uzPos,ndrange=(Nz,NumG))
end

# MomentumVectorInvariantCoriolisKernel! in the Cartesian frame. At the vertical node K of
# a layer the velocity is vK = uh + wK n (w of the lower/upper interface). The curl
#   J omega = sum_a D_a (J grad xi^a x v) + D_zeta (J grad zeta x v)
# is formed directly with the Cartesian metric rows, no rotation. The tendency -omega x vK
# goes into uh (projected after the assembly) and n . (-omega x vK) into w.
@kernel inbounds = true function MomentumVectorInvariantCoriolisCartKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(Glob),CoriolisFun,GridForm,Curvature,uzPos)

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  uCol = @localmem eltype(F) (N,N,3,ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+2)
  TZ = @localmem eltype(F) (N,N,3,2,ColumnTilesDim+2)
  YC = @localmem eltype(F) (N,N,3,2,2,ColumnTilesDim)

  if Iz <= Nz
    uCol[I,J,1,iz] = U[1,Iz,ind,2]
    uCol[I,J,2,iz] = U[1,Iz,ind,3]
    uCol[I,J,3,iz] = U[1,Iz,ind,uzPos]
    RhoCol[I,J,iz+1] = U[1,Iz,ind,1]
    wCol[I,J,iz+1] = U[1,Iz,ind,4]
    if iz == 1 && Iz == 1
      wCol[I,J,1] = BottomW(dXdxI,U[1,Iz,ind,2],U[1,Iz,ind,3],U[1,Iz,ind,uzPos],ID,IF)
      RhoCol[I,J,1] = U[1,Iz,ind,1]
    elseif iz == 1
      wCol[I,J,1] = U[1,Iz-1,ind,4]
      RhoCol[I,J,1] = U[1,Iz-1,ind,1]
    elseif iz == ColumnTilesDim && Iz < Nz  
      RhoCol[I,J,ColumnTilesDim+2] = U[1,Iz+1,ind,1]
    end   
  end 
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    n = SVector(FiniteElements.FrameNormal(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],GridForm))
    uh = SVector(uCol[I,J,1,iz],uCol[I,J,2,iz],uCol[I,J,3,iz])
    v1 = uh + wCol[I,J,iz] * n
    v2 = uh + wCol[I,J,iz+1] * n
#   J grad zeta x v at both vertical nodes, and the neighbours across the tile boundaries
    T1 = cross(MetricRowFull(dXdxI,3,1,ID,Iz,IF,n),v1)
    T2 = cross(MetricRowFull(dXdxI,3,2,ID,Iz,IF,n),v2)
    @unroll for d = 1 : 3
      TZ[I,J,d,1,iz+1] = T1[d]
      TZ[I,J,d,2,iz+1] = T2[d]
    end
    if iz == 1 && Iz > 1
      uhm1 = SVector(U[1,Iz-1,ind,2],U[1,Iz-1,ind,3],U[1,Iz-1,ind,uzPos])
      T = cross(MetricRowFull(dXdxI,3,2,ID,Iz-1,IF,n),uhm1 + wCol[I,J,1] * n)
      @unroll for d = 1 : 3
        TZ[I,J,d,2,1] = T[d]
      end
    end
    if iz == ColumnTilesDim && Iz < Nz
      uhp1 = SVector(U[1,Iz+1,ind,2],U[1,Iz+1,ind,3],U[1,Iz+1,ind,uzPos])
      T = cross(MetricRowFull(dXdxI,3,1,ID,Iz+1,IF,n),uhp1 + wCol[I,J,iz+1] * n)
      @unroll for d = 1 : 3
        TZ[I,J,d,1,iz+2] = T[d]
      end
    end
#   J grad xi^a x v, a = 1,2, for the horizontal derivatives
    @unroll for a = 1 : 2
      Y1 = cross(MetricRowFull(dXdxI,a,1,ID,Iz,IF,n),v1)
      Y2 = cross(MetricRowFull(dXdxI,a,2,ID,Iz,IF,n),v2)
      @unroll for d = 1 : 3
        YC[I,J,d,1,a,iz] = Y1[d]
        YC[I,J,d,2,a,iz] = Y2[d]
      end
    end
  end
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    n = SVector(FiniteElements.FrameNormal(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],GridForm))
    uh = SVector(uCol[I,J,1,iz],uCol[I,J,2,iz],uCol[I,J,3,iz])
    w1 = wCol[I,J,iz]
    w2 = wCol[I,J,iz+1]
    v1 = uh + w1 * n
    v2 = uh + w2 * n
    T1 = SVector(TZ[I,J,1,1,iz+1],TZ[I,J,2,1,iz+1],TZ[I,J,3,1,iz+1])
    T2 = SVector(TZ[I,J,1,2,iz+1],TZ[I,J,2,2,iz+1],TZ[I,J,3,2,iz+1])
    FH = SVector(eltype(F)(0),eltype(F)(0),eltype(F)(0))
    FW1 = eltype(F)(0)  
    FW2 = eltype(F)(0)  
#   jumps of J grad zeta x v at the interfaces
    if Iz > 1
      TB = SVector(TZ[I,J,1,2,iz],TZ[I,J,2,2,iz],TZ[I,J,3,2,iz])
      A = -cross(eltype(F)(0.5) * (T1 - TB),v1)
      FH += A
      FW1 += dot(n,A)
    end
    if Iz < Nz
      TT = SVector(TZ[I,J,1,1,iz+2],TZ[I,J,2,1,iz+2],TZ[I,J,3,1,iz+2])
      A = -cross(eltype(F)(0.5) * (TT - T2),v2)
      FH += A
      FW2 += dot(n,A)
    end
#   curl inside the layer at both vertical nodes
    c1 = D[I,1] * SVector(YC[1,J,1,1,1,iz],YC[1,J,2,1,1,iz],YC[1,J,3,1,1,iz]) +
      D[J,1] * SVector(YC[I,1,1,1,2,iz],YC[I,1,2,1,2,iz],YC[I,1,3,1,2,iz])
    c2 = D[I,1] * SVector(YC[1,J,1,2,1,iz],YC[1,J,2,2,1,iz],YC[1,J,3,2,1,iz]) +
      D[J,1] * SVector(YC[I,1,1,2,2,iz],YC[I,1,2,2,2,iz],YC[I,1,3,2,2,iz])
    for k = 2 : N
      c1 += D[I,k] * SVector(YC[k,J,1,1,1,iz],YC[k,J,2,1,1,iz],YC[k,J,3,1,1,iz]) +
        D[J,k] * SVector(YC[I,k,1,1,2,iz],YC[I,k,2,1,2,iz],YC[I,k,3,1,2,iz])
      c2 += D[I,k] * SVector(YC[k,J,1,2,1,iz],YC[k,J,2,2,1,iz],YC[k,J,3,2,1,iz]) +
        D[J,k] * SVector(YC[I,k,1,2,2,iz],YC[I,k,2,2,2,iz],YC[I,k,3,2,2,iz])
    end  
    Om = eltype(F)(0.5) * (T2 - T1)
    A1 = -cross(Om + c1,v1)
    A2 = -cross(Om + c2,v2)
    FH += A1 + A2
    FW1 += dot(n,A1)
    FW2 += dot(n,A2)

#   Shallow atmosphere: remove the tendencies uh*w/a, -|uh|^2/a of the shallow metric,
#   energy neutral in each half cell
    KinH = dot(uh,uh)
    FH -= Curvature * (JJ[ID,1,Iz,IF] * w1 + JJ[ID,2,Iz,IF] * w2) * uh
    FW1 += Curvature * JJ[ID,1,Iz,IF] * KinH
    FW2 += Curvature * JJ[ID,2,Iz,IF] * KinH

#   Coriolis in each half cell with its vertical velocity
    FhCor, FwCor = CoriolisFun(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],uh,w1)
    FH += FhCor * JJ[ID,1,Iz,IF]
    FW1 += FwCor * JJ[ID,1,Iz,IF]
    FhCor, FwCor = CoriolisFun(X[ID,2,1,Iz,IF],X[ID,2,2,Iz,IF],X[ID,2,3,Iz,IF],uh,w2)
    FH += FhCor * JJ[ID,2,Iz,IF]
    FW2 += FwCor * JJ[ID,2,Iz,IF]

    MC = M[Iz,ind,1] + M[Iz,ind,2]
    @atomic :monotonic F[1,Iz,ind,2] += FH[1] / MC
    @atomic :monotonic F[1,Iz,ind,3] += FH[2] / MC
    @atomic :monotonic F[1,Iz,ind,uzPos] += FH[3] / MC
    if Iz > 1  
      @atomic :monotonic F[1,Iz-1,ind,4] += RhoCol[I,J,iz+1] * FW1 /
        (RhoCol[I,J,iz] * M[Iz-1,ind,2] + RhoCol[I,J,iz+1] * M[Iz,ind,1])
    end
    if Iz < Nz
      @atomic :monotonic F[1,Iz,ind,4] += RhoCol[I,J,iz+1] * FW2 / 
        (RhoCol[I,J,iz+2] * M[Iz+1,ind,1] + RhoCol[I,J,iz+1] * M[Iz,ind,2])
    end  
  end
end

# sum of the metric row a over both vertical nodes of the layer
@inline MetricRowSum(dXdxI,a,ID,Iz,IF) = MetricRow(dXdxI,a,1,ID,Iz,IF) + MetricRow(dXdxI,a,2,ID,Iz,IF)

# Horizontal hyperviscosity in the Cartesian frame, first pass. The cache holds
#   1:3  n x grad(zeta), zeta the radial vorticity
#   4:6  grad(div uh)
#   7    div(grad theta)
# as Cartesian vectors. n . (J grad xi^a x uh) is formed at every node before the
# derivative, as Curl12 in the spherical frame.
@kernel inbounds = true function HyperViscCartKernel!(F,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),@Const(nR),uzPos) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  if Iz <= Nz
    n = SVector(nR[1,ind],nR[2,ind],nR[3,ind])
    uh = SVector(U[1,Iz,ind,2],U[1,Iz,ind,3],U[1,Iz,ind,uzPos])
    J1 = MetricRowSum(dXdxI,1,ID,Iz,IF)
    J2 = MetricRowSum(dXdxI,2,ID,Iz,IF)
    uCCol[I,J,iz] = dot(n,cross(J1,uh))
    vCCol[I,J,iz] = dot(n,cross(J2,uh))
    uDCol[I,J,iz] = dot(J1,uh)
    vDCol[I,J,iz] = dot(J2,uh)
    ThCol[I,J,iz] = U[1,Iz,ind,5] / U[1,Iz,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    Dxc = D[I,1] * ThCol[1,J,iz]
    Dyc = D[J,1] * ThCol[I,1,iz]
    Curl[I,J,iz] = D[I,1] * uCCol[1,J,iz] + D[J,1] * vCCol[I,1,iz] 
    Div[I,J,iz] = D[I,1] * uDCol[1,J,iz] + D[J,1] * vDCol[I,1,iz] 
    for k = 2 : N
      Dxc += D[I,k] * ThCol[k,J,iz]
      Dyc += D[J,k] * ThCol[I,k,iz] 
      Curl[I,J,iz] += D[I,k] * uCCol[k,J,iz] + D[J,k] * vCCol[I,k,iz] 
      Div[I,J,iz] += D[I,k] * uDCol[k,J,iz] + D[J,k] * vDCol[I,k,iz] 
    end
    Curl[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    Div[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    J1 = MetricRowSum(dXdxI,1,ID,Iz,IF)
    J2 = MetricRowSum(dXdxI,2,ID,Iz,IF)
    GradTh = (J1 * Dxc + J2 * Dyc) / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    ThCxCol[I,J,iz] = dot(J1,GradTh)
    ThCyCol[I,J,iz] = dot(J2,GradTh)
  end

  @synchronize 

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    DxCurl = DW[I,1] * Curl[1,J,iz]
    DyCurl = DW[J,1] * Curl[I,1,iz]
    DxDiv = DW[I,1] * Div[1,J,iz]
    DyDiv = DW[J,1] * Div[I,1,iz]
    DivTh = DW[I,1] * ThCxCol[1,J,iz] + DW[J,1] * ThCyCol[I,1,iz]
    for k = 2 : N
      DxCurl += DW[I,k] * Curl[k,J,iz]
      DyCurl += DW[J,k] * Curl[I,k,iz]
      DxDiv += DW[I,k] * Div[k,J,iz]
      DyDiv += DW[J,k] * Div[I,k,iz]
      DivTh += DW[I,k] * ThCxCol[k,J,iz] + DW[J,k] * ThCyCol[I,k,iz]
    end
    n = SVector(nR[1,ind],nR[2,ind],nR[3,ind])
    J1 = MetricRowSum(dXdxI,1,ID,Iz,IF)
    J2 = MetricRowSum(dXdxI,2,ID,Iz,IF)
    FC = cross(n,J1 * DxCurl + J2 * DyCurl)
    FD = J1 * DxDiv + J2 * DyDiv
    MC = M[Iz,ind,1] + M[Iz,ind,2]
    @unroll for d = 1 : 3
      @atomic :monotonic F[1,Iz,ind,d] += FC[d] / MC
      @atomic :monotonic F[1,Iz,ind,d+3] += FD[d] / MC
    end
    @atomic :monotonic F[1,Iz,ind,7] += DivTh / MC
  end
end

# second pass: curl of the curl part, divergence of the divergence part, tendency of uh
@kernel inbounds = true function HyperViscKoeffCartKernel!(F,@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffCurl,KoeffGrad,KoeffDiv,@Const(nR),uzPos)

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    n = SVector(nR[1,ind],nR[2,ind],nR[3,ind])
    C = SVector(Cache[1,Iz,ind,1],Cache[1,Iz,ind,2],Cache[1,Iz,ind,3])
    G = SVector(Cache[1,Iz,ind,4],Cache[1,Iz,ind,5],Cache[1,Iz,ind,6])
    J1 = MetricRowSum(dXdxI,1,ID,Iz,IF)
    J2 = MetricRowSum(dXdxI,2,ID,Iz,IF)
    uCCol[I,J,iz] = dot(n,cross(J1,C))
    vCCol[I,J,iz] = dot(n,cross(J2,C))
    uDCol[I,J,iz] = dot(J1,G)
    vDCol[I,J,iz] = dot(J2,G)
    ThCol[I,J,iz] = Cache[1,Iz,ind,7] 
  end
  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    Dxc = D[I,1] * ThCol[1,J,iz]
    Dyc = D[J,1] * ThCol[I,1,iz]
    Curl[I,J,iz] = D[I,1] * uCCol[1,J,iz] + D[J,1] * vCCol[I,1,iz]
    Div[I,J,iz] = D[I,1] * uDCol[1,J,iz] + D[J,1] * vDCol[I,1,iz]
    for k = 2 : N
      Dxc += D[I,k] * ThCol[k,J,iz]
      Dyc += D[J,k] * ThCol[I,k,iz]
      Curl[I,J,iz] += D[I,k] * uCCol[k,J,iz] + D[J,k] * vCCol[I,k,iz]
      Div[I,J,iz] += D[I,k] * uDCol[k,J,iz] + D[J,k] * vDCol[I,k,iz]
    end
    Curl[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    Div[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    J1 = MetricRowSum(dXdxI,1,ID,Iz,IF)
    J2 = MetricRowSum(dXdxI,2,ID,Iz,IF)
    GradTh = (J1 * Dxc + J2 * Dyc) / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    ThCxCol[I,J,iz] = U[1,Iz,ind,1] * dot(J1,GradTh)
    ThCyCol[I,J,iz] = U[1,Iz,ind,1] * dot(J2,GradTh)
  end

  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]
  if Iz <= Nz
    DxCurl = DW[I,1] * Curl[1,J,iz]
    DyCurl = DW[J,1] * Curl[I,1,iz]
    DxDiv = DW[I,1] * Div[1,J,iz]
    DyDiv = DW[J,1] * Div[I,1,iz]
    DivTh = DW[I,1] * ThCxCol[1,J,iz] + DW[J,1] * ThCyCol[I,1,iz]
    for k = 2 : N
      DxCurl += DW[I,k] * Curl[k,J,iz]
      DyCurl += DW[J,k] * Curl[I,k,iz]
      DxDiv += DW[I,k] * Div[k,J,iz]
      DyDiv += DW[J,k] * Div[I,k,iz]
      DivTh += DW[I,k] * ThCxCol[k,J,iz] + DW[J,k] * ThCyCol[I,k,iz]
    end
    n = SVector(nR[1,ind],nR[2,ind],nR[3,ind])
    J1 = MetricRowSum(dXdxI,1,ID,Iz,IF)
    J2 = MetricRowSum(dXdxI,2,ID,Iz,IF)
    FC = cross(n,J1 * DxCurl + J2 * DyCurl)
    FD = J1 * DxDiv + J2 * DyDiv
    FH = -(KoeffCurl * FC + KoeffGrad * FD) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,2] += FH[1]
    @atomic :monotonic F[1,Iz,ind,3] += FH[2]
    @atomic :monotonic F[1,Iz,ind,uzPos] += FH[3]
    @atomic :monotonic F[1,Iz,ind,5] += -KoeffDiv * DivTh / (M[Iz,ind,1] + M[Iz,ind,2])
  end
end
