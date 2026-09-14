import CGDycore:
  Thermodynamics, Sources, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DyCore, DGSEM
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using LinearAlgebra
using SparseArrays
import SparseArrays.ldiv!
using BandedMatrices
using FillArrays

function Permutation(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  for iz = 1 : nz
    for iv = 1 : 1
      for k = 1 : M
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
    for iv = 2 : 3
      for k = 2 : M - 1
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
  end
  ivw = 3
  ivTh = 2
  for iz = 1 : nz
      ii += 1
      p[ii] = 1 + (iz-1) * M  + (ivTh - 1) * N
      ii += 1
      p[ii] = 1 + (iz-1) * M  + (ivw - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivw - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
  end
  return p
end


function init_hydrostatic_perturbation(z; T0 = 300.0, z_center = 5000.0, width = 1000.0, dT = 1.0)
    p0 = 100000.0 # Pa
    g  = 9.81 # m/s^2
    R  = 287.05 # J/(kg·K)
    γ  = 1.4


    # 1. Hydrostatischer, isothermer Hintergrundzustand
    p_bg = p0 * exp(-g * z / (R * T0))
    ρ_bg = p_bg / (R * T0)
    θ_bg = T0 * (p0 / p_bg)^((γ - 1.0) / γ)

    # 2. Glatte Gauß-Störung auf die potentielle Temperatur
    Δθ = dT * exp(-((z - z_center)^2) / (2.0 * width^2))
    θ_total = θ_bg + Δθ
    ρ_total = ρ_bg

    return SVector{5, Float64}(ρ_total,0.0,0.0,0.0,ρ_total * θ_total)
end



function DScalarDMomAc(NZ,DG,cS)
  
  fac = 0.5
  M = DG.OrdPolyZ + 1
  N = NZ * M
  RowIndD = Int[]
  ColIndD = Int[]
  ValD = Float64[]
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  D = DG.DWZ
  for iZ = 1 : NZ
    for i = 1 : M
      for j = 1 : M
        push!(RowIndD,i+(iZ-1)*M)
        push!(ColIndD,j+(iZ-1)*M)
        push!(ValD,-D[i,j])
      end
    end  
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,fac/DG.wZ[1])
    end  
  end
  dSdM = sparse(RowInd, ColInd, Val,N,N)
  dSdG = sparse(RowIndD, ColIndD, ValD,N,N)

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/cS/DG.wZ[1])
    end  
  end
  dSdS = sparse(RowInd, ColInd, Val,N,N)
  return dSdS,dSdM,dSdG
end

function DMomDScalarAc(NZ,DG,cS)
  
  fac = 0.5
  M = DG.OrdPolyZ + 1
  N = NZ * M
  RowIndD = Int[]
  ColIndD = Int[]
  ValD = Float64[]
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  D = DG.DWZ
  for iZ = 1 : NZ
    for i = 1 : M
      for j = 1 : M
        push!(RowIndD,i+(iZ-1)*M)
        push!(ColIndD,j+(iZ-1)*M)
        push!(ValD,-D[i,j])
      end
    end  
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,fac/DG.wZ[1])
    end  
    if iZ == 1
      push!(RowInd,1)  
      push!(ColInd,1)  
      push!(Val,1/DG.wZ[1])
    end    
    if iZ == NZ
      push!(RowInd,M*NZ)  
      push!(ColInd,M*NZ)  
      push!(Val,-1/DG.wZ[1])
    end    
  end
  dMdG = sparse(RowIndD, ColIndD, ValD,N,N)
  dMdS = sparse(RowInd, ColInd, Val,N,N)
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac*cS/DG.wZ[1])
    end  
    if iZ == 1
      push!(RowInd,1)  
      push!(ColInd,1)  
      push!(Val,-cS/DG.wZ[1])
    end    
    if iZ == NZ
      push!(RowInd,M*NZ)  
      push!(ColInd,M*NZ)  
      push!(Val,-cS/DG.wZ[1])
    end    
  end  
  dMdM = sparse(RowInd, ColInd, Val,N,N)
  return dMdS,dMdM,dMdG
end

function DerivativeGeo(Geo,DG,dz)

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  DWS = DG.DWZ

  M = size(DWS,1)
  nz = size(dz,1)
  ND = size(dz,2)
  N = M * nz

  Gblk = zeros(M,M,nz,ND)

  @inbounds for ID in 1 : ND
    for iz in 1:nz
      inv2dz = 2 / dz[iz, ID]
      for i in 1:M
        acc = 0.0
        for k in 1:M
          dGeo = Geo[k,iz,ID] - Geo[i,iz,ID]
          Gblk[i,k,iz,ID] = inv2dz * 0.5 * DWS[i,k] * dGeo
          acc += DWS[i,k] * dGeo
          if ID == 1
            push!(RowInd,i+(iz-1)*M)
            push!(ColInd,k+(iz-1)*M)
            push!(Val,Gblk[i,k,iz,ID])
          end
        end
        if ID == 1
          push!(RowInd,i+(iz-1)*M)
          push!(ColInd,i+(iz-1)*M)
          push!(Val,inv2dz * 0.5 * acc)
        end
        Gblk[i,i,iz,ID] += inv2dz * 0.5 * acc
      end
    end
  end
  DGeo = sparse(RowInd,ColInd,Val,N,N)
  return Gblk, DGeo
end

function JacDGTNeu(U,Aux,DG,fac,dSdS,dSdM,dSdG,dMdS,dMdM,dMdG,dGeo,dz,Phys)
  D = DG.DWZ
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 5
  nz = size(U,2)
  M = size(U,1)
  oneM = ones(M)
  NF = size(dz,3)
  JacLU = Array{SparseArrays.UMFPACK.UmfpackLU}(undef,size(U,3))
    ID = 1
    @views dzCol = dz[:,ID]
    diagdz = spdiagm(2.0 ./ reshape(vec(oneM*dzCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Jac = [spzeros(N,N)  -diagdz* dSdS * diagm(dpdRhoTh) -diagdz * dSdM
           spzeros(N,N)  -diagdz*diagm(Th)*dSdS*diagm(dpdRhoTh)  -diagdz*diagm(Th)*dSdM
           spzeros(N,N)  -diagdz*dMdS*diagm(dpdRhoTh) -diagdz*dMdM]
    JacD = [sparse(fac*I,N,N)  spzeros(N,N) -diagdz * dSdG
           spzeros(N,N) sparse(fac*I,N,N)  -diagdz*dSdG*diagm(Th)
           dGeo -diagdz*dMdG*diagm(dpdRhoTh) sparse(fac*I,N,N)]
    JacLU[ID] = lu(Jac+JacD)
  return JacLU,Jac,JacD
end

function ExtendedJacobian(M,nz,dz,w,cS,Th,dpdRhoTh)
  N = 3 * M * nz

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]

  for iz = 1 : nz -1

#    Right part

    iRow = (iz - 1) * M + M 
    iCol = N + iz  
    ValLoc = 1.0 / (dz[iz]) 
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = iz * M + 1
    iCol = N + iz  
    ValLoc = -1.0 / (dz[iz+1]) 
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = M * nz + (iz - 1) * M + M 
    iCol = N + iz  
    ValLoc = (Th[(iz - 1) * M + M] + Th[iz * M + 1]) / (2 * dz[iz]) 
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = M * nz + iz * M + 1
    iCol = N + iz  
    ValLoc = -(Th[(iz - 1) * M + M] + Th[iz * M + 1]) / (2 * dz[iz+1]) 
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = 2 * M * nz + (iz - 1) * M + M 
    iCol = N + nz - 1 + iz + 1  
    ValLoc = -1.0 / (dz[iz]) 
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = 2 * M * nz + iz * M + 1
    iCol = N + nz - 1 + iz + 1  
    ValLoc = 1.0 / (dz[iz+1]) 
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

#   Low part

    iRow = N + iz  
    iCol = M * nz + (iz - 1) * M + M 
    ValLoc = -dpdRhoTh[(iz - 1) * M + M] / (cS * w)
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + iz  
    iCol = M * nz + iz * M + 1
    ValLoc = dpdRhoTh[iz * M + 1] / (cS * w)
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + iz  
    iCol = 2 * M * nz + (iz - 1) * M + M 
    ValLoc = -1.0 / w
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + iz  
    iCol = 2 * M * nz + iz * M + 1
    ValLoc = -1.0 / w
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + nz - 1 + iz + 1 
    iCol = M * nz + (iz - 1) * M + M 
    ValLoc = dpdRhoTh[(iz - 1) * M + M] / (w)
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + nz - 1 + iz + 1 
    iCol = M * nz + iz * M + 1
    ValLoc = dpdRhoTh[iz * M + 1] / (w)
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + nz - 1 + iz + 1  
    iCol = 2 * M * nz + (iz - 1) * M + M 
    ValLoc = cS / w
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)

    iRow = N + nz - 1 + iz + 1   
    iCol = 2 * M * nz + iz * M + 1
    ValLoc = -cS / w
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)
  end  

# Right part
  iRow = 2 * M * nz + 1
  iCol = N + nz  
  ValLoc = 2.0 / (dz[1]) 
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,ValLoc)

  iRow = 3 * M * nz 
  iCol = N + 2 * nz 
  ValLoc = -2.0 / (dz[nz]) 
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,ValLoc)


# Low part

  iRow = N + nz   
  iCol = 2 * M * nz + 1
  ValLoc = -1.0 * cS / w
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,ValLoc)

  iRow = N + 2 * nz 
  iCol = 3 * M * nz 
  ValLoc = 1.0 * cS / w
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,ValLoc)

  iRow = N + nz 
  iCol = M * nz + 1
  ValLoc = 1.0 * dpdRhoTh[1] / w
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,ValLoc)

  iRow = N + 2 * nz
  iCol = 2 * M * nz 
  ValLoc = 1.0 * dpdRhoTh[end] / w
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,ValLoc)
  for iz = 1 : 2 * nz 
    iRow = N + iz  
    iCol = N + iz
    ValLoc = 1.0
    push!(RowInd,iRow)
    push!(ColInd,iCol)
    push!(Val,ValLoc)
  end  

  JacExtend = sparse(RowInd,ColInd,Val,N + 2 * nz,N + 2 * nz)
end

backend = CPU()
FTB = Float64
Parallel = true

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

Problem = "WarmBubble2DXCart"
Param = Examples.Parameters(FTB,Problem)
# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

# ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Grid
Boundary = Grids.Boundary()
Boundary.WE = "Period"
Boundary.SN = "Period"
Boundary.BT = "FreeSlip"
nx = 2
ny = 2
Lx = 20000.0
Ly = 2000.0
H = 10000.0
x0 = 0.0
y0 = 0.0
nz = 10
OrdPoly = 4
OrdPolyZ = 6
M = OrdPolyZ + 1
N = M * nz
Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)
Grid.AdaptGrid = Grids.AdaptGrid(FTB,"Sleve",FTB(H))

OrdPrint = 4
OrdPrintZ = 4
DGMethod = ""
TopoS = ""
Topography=(TopoS=TopoS,
              )
TopoProfile = Examples.Flat()()
(DG, Metric, Exchange, Global) = DyCore.InitCartDG(backend,FTB,OrdPoly,OrdPolyZ,DGMethod,
    OrdPrint,OrdPrintZ,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)
#DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)


fac = 0.5
U = ones(M,nz,DG.NumI,5)
Aux = rand(M,nz,DG.NumI,4)
@. Aux += 50.0
H = 10000.0
dzLoc = H / nz
dz = ones(nz,DG.NumI) * dzLoc
dz = ones(nz,DG.NumI) * dzLoc
dz[10] = 1450.0
dz[9] = 1350.0
dz[8] = 1250.0
dz[7] = 1150.0
dz[6] = 1050.0
dz[5] = 950.0
dz[4] = 850.0
dz[3] = 750.0
dz[2] = 650.0
dz[1] = 550.0
Metric.dz = dz

VelForm = Examples.VelocityC()
time = 0.0
for iz = 1 : nz
  z0 = (iz - 1) * dzLoc  
  z1 = iz * dzLoc  
  for k = 1 : M
    zLoc = 0.5 * ((1.0 - DG.xwZ[k]) * z0 + (1.0 + DG.xwZ[k]) * z1)  
    xS = SVector{3}(0.0,0.0,zLoc)
    RhoP,_,_,_,RhoThP= init_hydrostatic_perturbation(zLoc)
    @views @. U[k,iz,:,1] = RhoP*(1+0.1*rand())
    @views @. U[k,iz,:,5] = RhoThP + 2.0 * rand()
  end
end

Model.GPAuxPos = 2
Model.GeoPotential = Sources.GeoPotentialDeep()(Model.GPAuxPos,Grid.Form)
DGSEM.GeoPot(Aux,DG,Metric,Exchange,Global)
Geo = Aux[:,:,:,Model.GPAuxPos]

Gblk, dGeo = DerivativeGeo(Geo,DG,dz)
dSdS,dSdM,dSdG = DScalarDMomAc(nz,DG,Param.cS)
dMdS,dMdM,dMdG = DMomDScalarAc(nz,DG,Param.cS)

fac = 0.5
invfac=1.0/fac
JacLU, Jac, JacD =JacDGTNeu(U,Aux,DG,invfac,dSdS,dSdM,dSdG,dMdS,dMdM,dMdG,dGeo,dz,Phys)
dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
  (Phys.Rd * U[:,:,1,5] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),M * nz)
Th = reshape(U[:,:,1,5] ./ U[:,:,1,1], M * nz)

JacExt = ExtendedJacobian(M,nz,dz,DG.wZ[1],Param.cS,Th,dpdRhoTh)

n11 = 3 * M * nz
JacExt11 = JacExt[1:n11,1:n11]
JacExt12 = JacExt[1:n11,n11+1:end]
JacExt21 = JacExt[n11+1:end,1:n11]
JacExt22 = JacExt[n11+1:end,n11+1:end]

JacExt11F = collect(JacExt11)
JacExt12F = collect(JacExt12)
JacExt21F = collect(JacExt21)
JacExt22F = collect(JacExt22)

SchurJacExtF = JacExt11F - JacExt12F * (JacExt22F \ JacExt21F)
SchurJacExt = sparse(SchurJacExtF)

JacExtD = deepcopy(JacExt)
JacExtD[1:n11,1:n11] += JacD

JacExtD11 = JacExtD[1:n11,1:n11]
@show JacExtD11[1,1]
JacExtD12 = JacExtD[1:n11,n11+1:end]
JacExtD21 = JacExtD[n11+1:end,1:n11]
JacExtD22 = JacExtD[n11+1:end,n11+1:end]

JacExtD11F = collect(JacExtD11)
JacExtD12F = collect(JacExtD12)
JacExtD21F = collect(JacExtD21)
JacExtD22F = collect(JacExtD22)

SchurJacExtDF = JacExtD22F - JacExtD21F * (JacExtD11F \ JacExtD12F)
SchurJacExtD = sparse(SchurJacExtDF)

pS = zeros(Int,2*nz)
pS[1] = nz 
ii = 1
for i = 1 : nz - 1
  pS[ii+1] = i
  pS[ii+2] = nz  + i
  global ii += 2
end
pS[2*nz] = 2*nz
SchurJacExtDPF = SchurJacExtDF[pS,pS]
SchurJacExtDP = SchurJacExtD[pS,pS]
SchurJacExtDPB = BandedMatrix(SchurJacExtDP)

# Linear algebra test

b = ones(size(U))
@. b[:,:,:,5] *= 2
@. b[:,:,:,4] *= 3
bSplit = similar(b)
@. bSplit = b

bF = zeros(M*nz*3)
bHDGF = zeros(M*nz*3)
bFExt = zeros(M*nz*3+2*nz)
ibF = 0
@inbounds for iv in [1,5,4]
  @inbounds for iz = 1 : nz
    @inbounds for i = 1 : M
      global ibF += 1
      bF[ibF] = b[i,iz,1,iv]
      bFExt[ibF] = b[i,iz,1,iv]
    end
  end
end



p1=collect([collect(1:7);collect(71:77);collect(141:147)])
p2=p1.+M
p3=p1.+2*M
x = (Jac + JacD) \ bF
xExt = JacExtD \ bFExt

y1 = JacExtD11F \ bF
y1p1 = y1[p1]
@show y1p1[8:14]
@show y1p1[15:21]
y2Vor = JacExtD21F * y1
y2 = -SchurJacExtDF\(JacExtD21F * y1)
y1 = JacExtD11F \ (bF - JacExtD12F * y2)

JacHDG = DGSEM.JacHDGVert(backend,FTB,M,nz,DG)
DGSEM.Jac!(U,fac,DG,Metric,Phys,Aux,JacHDG,Global,VelForm)
DGSEM.Solve!(JacHDG,b)
JacSplit = DGSEM.JacSplitDGVert(backend,FTB,M,nz,DG)
DGSEM.Jac!(U,fac,DG,Metric,Phys,Aux,JacSplit,Global,VelForm)
DGSEM.Solve!(JacSplit,bSplit)
ldiv!(JacLU[1],bF)
stop

pE=p1.+(nz-1)*M
@show p1
@show p2
pp1=[p1;collect(211:230)]
pp2=[p2;collect(211:230)]

JacExtD11F_1 = JacExtD11F[p1,p1]
JacExtD11F_2 = JacExtD11F[p2,p2]
JacExtD11F_3 = JacExtD11F[p3,p3]
JacExtD11F_E = JacExtD11F[pE,pE]
r1 = JacExtD12F[p1,:]
r2 = JacExtD12F[p2,:]
r3 = JacExtD12F[p3,:]
r1_1 = r1[1:7,pS[1]]
r1_2 = r1[8:14,pS[1]]
r1_3 = r1[15:21,pS[1]]
@show r1_1
@show r1_2
@show r1_3
rE = JacExtD12F[pE,:]
c1 = JacExtD21F[:,p1]
c2 = JacExtD21F[:,p2]
c3 = JacExtD21F[:,p3]
c1_1 = c1[pS[1],1:7]
c1_2 = c1[pS[1],8:14]
c1_3 = c1[pS[1],15:21]
@show c1_1
@show c1_2
@show c1_3
cE = JacExtD21F[:,pE]
rr1 = (JacExtD11F_1 \ r1)
rr1_1 = rr1[1:7,pS[1]]
rr1_2 = rr1[8:14,pS[1]]
rr1_3 = rr1[15:21,pS[1]]
@show rr1_1
@show rr1_2
@show rr1_3
S1 =  JacExtD22F - c1 * (JacExtD11F_1 \ r1)
S1P = S1[pS,pS]
S2 =  JacExtD22F - c2 * (JacExtD11F_2 \ r2)
S2P = S2[pS,pS]
luJacExtD11F_3 = lu(JacExtD11F_3,NoPivot())
S3 =  JacExtD22F - c3 * (JacExtD11F_3 \ r3)
S3P = S3[pS,pS]
r3I = deepcopy(r3)
r3II = -(JacExtD11F_3[15:21,1:7]*r3_1 + JacExtD11F_3[15:21,8:14]*r3_2) / fac
ldiv!(luJacExtD11F_3,r3I)
@show r3I[:,pS[4]]
S3lu = JacExtD22F - c3 * r3I
S3luP = S3lu[pS,pS]
SE =  JacExtD22F - cE * (JacExtD11F_E \ rE)
SEP = SE[pS,pS]
stop
