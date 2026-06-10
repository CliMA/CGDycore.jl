import CGDycore:
  Parameters, Thermodynamics, Examples, Sources, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGSEM, CGSEM, DyCore, IMEXRosenbrock

using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll
using MPI

#=
abstract type RiemannSolver end

Base.@kwdef struct RiemannLMARSLinFast <: RiemannSolver end
function (::RiemannLMARSLinFast)(Phys,RhoPos,uPos,vPos,wPos,RhoThPos,dpdRhoThPos,ThPos)
  @inline function RiemannByLMARSSemi!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    pLL = AuxL[dpdRhoThPos] * VLL[RhoThPos]
    pRR = AuxR[dpdRhoThPos] * VRR[RhoThPos]
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3)
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3)
    pM = FT(0.5) * ((pLL + pRR) - Phys.cS * (vRR - vLL))
    vM = FT(0.5) * ((vRR + vLL) - Phys.invcS * (pRR - pLL))
    F[uPos] = n1 * pM
    F[vPos] = n2 * pM
    F[wPos] = n3 * pM
    F[RhoPos] = vM
    F[RhoThPos] = FT(0.5) * F[RhoPos] * (AuxL[ThPos] + AuxR[ThPos])
  end
  return RiemannByLMARSSemi!
end

@kernel inbounds = true function RiemannNonLinV3Kernel!(RiemannSolver!,F,@Const(U),
  @Const(Aux),@Const(NV),@Const(VolSurfV),
  ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  Iz,ind = @index(Global, NTuple)


  Nz = @uniform @ndrange()[1]
  NumI = @uniform @ndrange()[2]

  VLL = @private eltype(F) (NUMV,)
  VRR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)
  FLoc = @private eltype(F) (NUMV,)

  if ind <= NumI
    n1 = NV[Iz,ind,1]
    n2 = NV[Iz,ind,2]
    n3 = NV[Iz,ind,3]
    if Iz == 1
      @unroll for iAux = 1 : NAUX
        AuxL[iAux] = Aux[1,Iz,ind,iAux]
        AuxR[iAux] = AuxL[iAux]
      end
      @unroll for iv = 1 : NUMV
        VLL[iv] = U[1,Iz,ind,iv]
        VRR[iv] = VLL[iv]
      end
      t = eltype(F)(2) * (n1 * VLL[2] +
        n2 * VLL[3] +
        n3 * VLL[4])
      VLL[2] -= n1 * t
      VLL[3] -= n2 * t
      VLL[4] -= n3 * t
      RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]  
      @unroll for iv = 1 : NUMV  
        F[1,Iz,ind,iv] += FLoc[iv] * Surf
      end  
    elseif Iz == Nz 
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[end,Iz-1,ind,iAux]
        AuxR[iAux] = AuxL[iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[end,Iz-1,ind,iv]
        VRR[iv] = VLL[iv]
      end  
      t = eltype(F)(2) * (n1 * VRR[2] +
        n2 * VRR[3] +
        n3 * VRR[4])
      VRR[2] -= n1 * t
      VRR[3] -= n2 * t
      VRR[4] -= n3 * t
      RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]   
      @unroll for iv = 1 : NUMV  
        F[end,Iz-1,ind,iv] -= FLoc[iv] * Surf
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[end,Iz-1,ind,iAux]
        AuxR[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[end,Iz-1,ind,iv]
        VRR[iv] = U[1,Iz,ind,iv]
      end  
      RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]  
      @unroll for iv = 1 : NUMV  
        FLoc[iv] *= Surf
        F[end,Iz-1,ind,iv] -= FLoc[iv]
        F[1,Iz,ind,iv] += FLoc[iv]
      end  
    end    
  end  
end
=#

function DLagrange(x,xP,iP)
  Df=eltype(x)(0)
  f=eltype(x)(1)
  for i = 1 : size(xP,1)
    if i ≠ iP
      Df = Df * (x-xP[i]) / (xP[iP] - xP[i]) + f / (xP[iP]-xP[i])
      f = f * (x - xP[i]) / (xP[iP] - xP[i])
    end
  end
  return Df
end

###### Main




  backend = CPU()
  FTB = Float64
  MPI.Init()
  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = DyCore.ParallelComStruct()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  Phys = DyCore.PhysParameters{FTB}()
  Model = DyCore.ModelStruct{FTB}()
  OrdPoly = 4
  OrdPolyZ = 4
  N = OrdPoly + 1
  M = OrdPolyZ + 1
  nz = 13
  nPanel = 5
  RefineLevel = 3
  ns = 10
  nLon = 10
  nLat = 10
  LatB = 80
  RadEarth = 1.0
  GridType = "CubedSphere"
  Decomp="EqualArea"
  Discretization="DG"
  DGMethod = "Kubatko2LGL"
  AdaptGridType = "Sleve"
  H = 1000.0
  OrdPrint = 1
  OrdPrintZ = 1
  TopoS = ""
  TopoProfile = Examples.Flat()()
  Problem = "BaroWaveDrySphere"

  Param = Examples.Parameters(FTB,Problem)

  Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,
    GridType,Decomp,RadEarth,Model,ParallelCom;Discretization=Discretization,ChangeOrient=2)

  Grid.AdaptGrid = Grids.AdaptGrid(FTB,AdaptGridType,FTB(H))
  Topography = (TopoS=TopoS,H=H,Rad=RadEarth)

  (DG, Metric, Exchange, Global) = DyCore.InitSphereDG(backend,FTB,OrdPoly,OrdPolyZ,DGMethod,
    OrdPrint,OrdPrintZ,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5


  NUMV = 5
  nz = Grid.nz
  NF = size(DG.Glob,2)
  F = zeros(M,nz,DG.NumI,NUMV)
  FRV = zeros(M,nz,DG.NumI,NUMV)
  FRH = zeros(M,nz,DG.NumI,NUMV)
  U = rand(M,nz,DG.NumI,NUMV)
  NAUX = 4
  Aux = rand(M,nz,DG.NumI,NAUX)
  @. Aux[:,:,:,2] += 10 
  @. Aux[:,:,:,3] += 400
  @. Aux[:,:,:,4] += 300

  JacVH = DGSEM.JacFluxVolumeSparseH(DG,Metric,Aux,Phys)
  JacVV = DGSEM.JacFluxVolumeSparseV(DG,Metric,Aux,Phys)
  JacRH = DGSEM.JacRiemannSparseH(DG,Metric,Grid,Aux,Phys)
  JacRV = DGSEM.JacRiemannSparseV(DG,Metric,Aux,Phys)
  JacGeoH = DGSEM.JacGravitySparseH(DG,Metric,Aux,Phys)
  JacGeoV = DGSEM.JacGravitySparseV(DG,Metric,Aux,Phys)

  U1 = reshape(U,M*nz*DG.NumI*NUMV)
  FVH = JacVH * U1
  FVV = JacVV * U1
  FRH = JacRH * U1
  FRV = JacRV * U1
  FGeoH = JacGeoH * U1
  FGeoV = JacGeoV * U1

  FM = FVH + FVV + FRH + FRV + FGeoH + FGeoV
  F1 = reshape(FM,M,nz,N,N,NF,NUMV)

  RiemannSolverFast = DGSEM.RiemannLMARSLinFast()(Phys,Phys,RhoPos,uPos,vPos,
    wPos,RhoThPos,3,4)

  FluxAverageFast = DGSEM.KennedyGruberGravLinFast1()(RhoPos,uPos,vPos,wPos,
    RhoThPos,3,4,2,Grid.Type)

  group = (nz+1,1)
  ndrange = (nz+1,DG.NumI)
  KRiemannNonLinV3Kernel! = DGSEM.RiemannNonLinV3Kernel!(backend,group)
  KRiemannNonLinV3Kernel!(RiemannSolverFast,F,U,Aux,Metric.NV,
    Metric.VolSurfV,Val(NUMV),Val(NAUX);ndrange=ndrange)

  DoFE = DG.DoFE
  NE = Grid.NumEdges
  group = (M,nz,1,1)
  ndrange = (M,nz,DoFE,NE)
  KRiemannNonLinH3Kernel! = DGSEM.RiemannNonLinH3Kernel!(backend,group)
  KRiemannNonLinH3Kernel!(RiemannSolverFast,F,U,Aux,DG.GlobE,Grid.EF,Grid.FE,Metric.NH,
    Metric.VolSurfH,DG.wF,Grid.NumFaces,Val(NUMV),Val(NAUX);ndrange=ndrange)

  group = (N,N,M,1,1)
  ndrange = (N,N,M,nz,NF)
  KFluxSplitVolumeNonLinHKernel! = DGSEM.FluxSplitVolumeNonLinHVQuadKernel!(backend,group)
  @views KFluxSplitVolumeNonLinHKernel!(FluxAverageFast,F,U,Aux,Metric.dXdxI,DG.DVT,DG.DVZT,DG.Glob,
    Val(N), Val(M), Val(NUMV), Val(NAUX);ndrange=ndrange)
  F2 = reshape(F,M,nz,N,N,NF,NUMV)

  permu = Int[]
  for ind = 1 : DG.NumI
      for iz = 1 : nz  
        for k = 1 : M
    for iv = 1 : NUMV  
          ip = k + (iz - 1) * M + (ind - 1) * M *nz + (iv - 1) * M * nz * DG.NumI  
          push!(permu,ip)
        end
      end
    end
  end  
  JacV = JacVV + JacRV + JacGeoV
  JacVP = JacV[permu,permu]
  nn = M*nz*NUMV
  B = JacVP[1:nn,1:nn]
  aa = 3;


