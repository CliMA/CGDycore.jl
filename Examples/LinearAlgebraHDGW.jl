import CGDycore:
  Parameters as P, Thermodynamics, Sources, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DyCore, DGSEM
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
OrdPolyZ = 4
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
    @views @. U[k,iz,:,1] = RhoP #*(1+0.1*rand())
    @views @. U[k,iz,:,5] = RhoThP  #+ 2.0 * rand()
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
pp = Permutation(M,nz)
pBl = zeros(Int,3*M*nz)
pB = zeros(Int,3*M*nz)
pR = zeros(Int,2*M*nz)
pRL = zeros(Int,2*M*nz+nz-1)
ii1 = 0
iB = 0
iB1 = 0
ii2 = M * nz
ii3 = 2 * M * nz
for iz = 1 : nz
  for k = 1 : M
    global ii1 += 1
    pBl[ii1] = k + (iz - 1) * 3 * M
    global ii2 += 1
    pBl[ii2] = k + M + (iz - 1) * 3 * M
    global ii3 += 1
    pBl[ii3] = k + 2 * M + (iz - 1) * 3 * M
  end
  for k = 1 : M
    global iB += 1  
    pB[iB] = k + (iz - 1) * M   
  end
  for k = 1 : M
    global iB += 1  
    global iB1 += 1  
    pB[iB] = k + (iz - 1) * M + M * nz  
    pR[iB1] = iB
  end
  for k = 1 : M
    global iB += 1  
    global iB1 += 1  
    pB[iB] = k + (iz - 1) * M + 2 * M * nz  
    pR[iB1] = iB
  end
end
pRL[1:2*M*nz] .= pR
for iz = 1 : nz -1
  global iB1 += 1
  pRL[iB1] = 3 * M * nz + iz
end  
JacP = Jac[pB,pB]
JacDP = JacD[pB,pB]
JacAllP = JacP + JacDP

# Linear algebra test

b = ones(size(U))
@. b[:,:,:,5] *= 2
@. b[:,:,:,4] *= 3
bSplit = similar(b)
@. bSplit = b

bF = zeros(M*nz*3)
bHDGF = zeros(M*nz*3)
bFExt = zeros(M*nz*3+nz-1)
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


JacHDG = DGSEM.JacHDGVert(backend,FTB,M,nz,DG)
DGSEM.Jac!(U,fac,DG,Metric,Phys,Aux,JacHDG,Global,VelForm)
DGSEM.Solve!(JacHDG,b)


function JacSparseHDG(Jac)
  M = JacHDG.M
  nz = JacHDG.nz
  ID = 1
  fac = Jac.fac
  @show "HDG",fac


  # Pre-allocate array of sparse matrices for each iz slice
  sparse_blocks11 = Vector{SparseMatrixCSC{Float64, Int}}(undef, nz)
  sparse_blocks12 = Vector{SparseMatrixCSC{Float64, Int}}(undef, nz)
  sparse_blocks21 = Vector{SparseMatrixCSC{Float64, Int}}(undef, nz)
  sparse_blocks22 = Vector{SparseMatrixCSC{Float64, Int}}(undef, nz-1)

  b1 = zeros(M,nz-1)
  b2 = zeros(M,nz-1)
  b3 = zeros(M,nz-1)
  c1 = zeros(nz-1,M)
  c2 = zeros(nz-1,M)
  c3 = zeros(nz-1,M)
  for iz in 1:nz
    # Extract 2D slices for current iz
    a13 = Jac.A13[:, :, iz, ID]
    a23 = Jac.A23[:, :, iz, ID]
    a31 = Jac.A31[:, :, iz, ID]
    a32 = Jac.A32[:, :, iz, ID]
    a33 = I(M)/fac
    a33[1,1] += Jac.A33[iz,ID]
    a33[M,M] += Jac.A33[iz,ID]
    # Assemble full 2D block slice
    block11 = [I(M)/fac        zeros(M, M) a13;
               zeros(M, M) I(M)/fac        a23;
               a31         a32             a33]
    # Convert slice directly to sparse
    sparse_blocks11[iz] = sparse(block11)
  end             

  for iz = 1 : nz
    @. b1 = 0.0    
    @. b2 = 0.0    
    @. b3 = 0.0    
    @. c2 = 0.0    
    @. c3 = 0.0    
    if iz < nz
      b1[M,iz] = Jac.B1[2,iz,ID]  
      b2[M,iz] = Jac.B2[2,iz,ID]  
      b3[M,iz] = Jac.B3[2,iz,ID]  
      c2[iz,M] = Jac.C2[2,iz,ID]  
      c3[iz,M] = Jac.C3[2,iz,ID]  
      d = 2.0 * P.cS
      sparse_blocks22[iz] = sparse([d])
    end    
    if iz > 1
      b1[1,iz-1] = Jac.B1[1,iz,ID]  
      b2[1,iz-1] = Jac.B2[1,iz,ID]  
      b3[1,iz-1] = Jac.B3[1,iz,ID]  
      c2[iz-1,1] = Jac.C2[1,iz,ID]  
      c3[iz-1,1] = Jac.C3[1,iz,ID]  
    end    
    block12 = [b1
               b2
               b3]
    block21 = [c1 c2 c3]
    sparse_blocks12[iz] = sparse(block12)
    sparse_blocks21[iz] = sparse(block21)
  end

  b11S = blockdiag(sparse_blocks11...)
  b12S = vcat(sparse_blocks12...)
  b21S = hcat(sparse_blocks21...)
  b22S = blockdiag(sparse_blocks22...)

  # Efficiently construct sparse block diagonal matrix
  return [b11S b12S
          b21S b22S]
end
JacHDGS = JacSparseHDG(JacHDG)
JacHDGSR = JacHDGS[pR,pR]

JacSplit = DGSEM.JacSplitDGVert(backend,FTB,M,nz,DG)
DGSEM.Jac!(U,fac,DG,Metric,Phys,Aux,JacSplit,Global,VelForm)
@show "vor Solve JacSplit"
DGSEM.Solve!(JacSplit,bSplit)

function JacSparseSplit(Jac)
  M = JacHDG.M
  nz = JacHDG.nz
  ID = 1
  fac = Jac.fac
  @show "Split",fac


  # Pre-allocate array of sparse matrices for each iz slice
  sparse_blocks = Vector{SparseMatrixCSC{Float64, Int}}(undef, nz)
  for iz = 1 : nz
    a13 = [Jac.B1_23[:,1,iz,ID] Jac.A13[:, :, iz, ID]  Jac.B1_23[:,2,iz,ID]]
    a23 = [Jac.DD[ID,1,2,iz]    Jac.C14_3[1:1,:,iz,ID] Jac.DD[ID,1,3,iz]
           Jac.B2_23[:,1,iz,ID] Jac.A23[:, :, iz, ID]  Jac.B2_23[:,2,iz,ID]
           Jac.DD[ID,4,2,iz]    Jac.C14_3[2:2,:,iz,ID] Jac.DD[ID,4,3,iz]]
    a31 = [Jac.C23_1[1:1,:,iz,ID]
           Jac.A31[:,:,iz,ID]
           Jac.C23_1[2:2,:,iz,ID]]
    a32 = [Jac.DD[ID,2,1,iz] Jac.C23_2[1:1,:,iz,ID] Jac.DD[ID,2,4,iz]
           Jac.B3_14[:,1:1,iz,ID] Jac.A32[:,:,iz,ID] Jac.B3_14[:,2:2,iz,ID]
           Jac.DD[ID,3,1,iz]  Jac.C23_2[2:2,:,iz,ID] Jac.DD[ID,3,4,iz] ]
    # Assemble full 2D block slice
    block = [I(M)*fac    zeros(M, M)  a13;
               zeros(M, M) I(M)*fac     a23;
               a31         a32          I(M)*fac]
    block[M,2*M] = Jac.B1_4[iz,ID]           
    block[M,3*M] = Jac.B1p_12[2,iz,ID]
    block[2*M,2*M] = Jac.DD[ID,4,4,iz]
#   if iz < nz
#     block[2*M,3*M] = Jac.DU[ID,4,2,iz]
#   end  
    # Convert slice directly to sparse
    sparse_blocks[iz] = sparse(block)
  end

  bS = blockdiag(sparse_blocks...)
end
JacSplitS = JacSparseSplit(JacSplit)



n11 = 3 * M * nz
JacHDGF11 = collect(JacHDGS[1:n11,1:n11])
JacHDGF12 = collect(JacHDGS[1:n11,n11+1:end])
JacHDGF21 = collect(JacHDGS[n11+1:end,1:n11])
JacHDGF22 = collect(JacHDGS[n11+1:end,n11+1:end])
JacHDGFSchur = JacHDGF11 - JacHDGF12 * (JacHDGF22 \ JacHDGF21)
JacHDGSSchur = sparse(JacHDGFSchur)
JacHDGFSchurR = JacHDGF22 - JacHDGF21 * (JacHDGF11 \ JacHDGF12)
stop
