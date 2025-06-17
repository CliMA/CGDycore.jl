import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGVertical, GPU, DyCore, DGSEM
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

mutable struct Interior 
  M::Int
  nz::Int
  fac::Float64
  FacGrav::Float64
  A13::Array{Float64, 4}
  A23::Array{Float64, 4}
  A31::Array{Float64, 4}
  A32::Array{Float64, 4}
  B13::Array{Float64, 4}
  B23::Array{Float64, 4}
  B32::Array{Float64, 4}
  C13::Array{Float64, 4}
  C23::Array{Float64, 4}
  C32::Array{Float64, 4}
  SA::Array{Float64, 4}
end  

function Interior(M,nz,NumG)
  M2 = M - 2
  fac = 0
  FacGrav = 0
  A13 = zeros(M2,M2,nz,NumG)
  A23 = zeros(M2,M2,nz,NumG)
  A31 = zeros(M2,M2,nz,NumG)
  A32 = zeros(M2,M2,nz,NumG)
  B13 = zeros(M2,2,nz,NumG)
  B23 = zeros(M2,2,nz,NumG)
  B32 = zeros(M2,2,nz,NumG)
  C13 = zeros(2,M2,nz,NumG)
  C23 = zeros(2,M2,nz,NumG)
  C32 = zeros(2,M2,nz,NumG)
  SA = zeros(M2,M2,nz,NumG)

  return Interior(
    M,
    nz,
    fac,
    FacGrav,
    A13,
    A23,
    A31,
    A32,
    B13,
    B23,
    B32,
    C13,
    C23,
    C32,
    SA,
  )
end  

function FillInteriorP1!(AI,JacP11,fac,Phys)
  M2 = AI.M - 2
  M23= M2 * 3
  nz = AI.nz
  AI.fac = fac
  AI.FacGrav = 0.5 * Phys.Grav

  for iz = 1 : nz
    A = JacP11[1+(iz-1)*M23:iz*M23,1+(iz-1)*M23:iz*M23]
    @views @. AI.A13[:,:,iz,1] = A[1:M2,2*M2+1:3*M2]
    @views @. AI.A23[:,:,iz,1] = A[M2+1:2*M2,2*M2+1:3*M2]
    @views @. AI.A31[:,:,iz,1] = A[2*M2+1:3*M2,1:M2]
    @views @. AI.A32[:,:,iz,1] = A[2*M2+1:3*M2,M2+1:2*M2]
    @views AI.SA[:,:,iz,1] = fac * I - (0.5 * Phys.Grav / fac) * AI.A13[:,:,iz,1] - 
      (1.0 / fac) *AI.A23[:,:,iz,1] * AI.A32[:,:,iz,1]
    B = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-1)*6:iz*6]
    @views @. AI.B13[:,:,iz,1] = B[1:M2,5:6]
    @views @. AI.B23[:,:,iz,1] = B[M2+1:2*M2,5:6]
    @views @. AI.B32[:,:,iz,1] = B[2*M2+1:3*M2,3:4]
    C = JacP21[1+(iz-1)*6:iz*6,1+(iz-1)*M23:iz*M23]
    @views @. AI.C13[:,:,iz,1] = C[1:2,2*M2+1:3*M2]
    @views @. AI.C23[:,:,iz,1] = C[3:4,2*M2+1:3*M2]
    @views @. AI.C32[:,:,iz,1] = C[5:6,M2+1:2*M2]
  end
end

function FillInteriorP2!(AI,JacP11,fac,Phys)
  M2 = AI.M - 2
  M23= M2 * 3
  nz = AI.nz
  AI.fac = fac
  AI.FacGrav = 0.5 * Phys.Grav

  for iz = 1 : nz
    A = JacP11[1+(iz-1)*M23:iz*M23,1+(iz-1)*M23:iz*M23]
    @views @. AI.A13[:,:,iz,1] = A[1:M2,2*M2+1:3*M2]
    @views @. AI.A23[:,:,iz,1] = A[M2+1:2*M2,2*M2+1:3*M2]
    @views @. AI.A31[:,:,iz,1] = A[2*M2+1:3*M2,1:M2]
    @views @. AI.A32[:,:,iz,1] = A[2*M2+1:3*M2,M2+1:2*M2]
    @views AI.SA[:,:,iz,1] = fac * I - (0.5 * Phys.Grav / fac) * AI.A13[:,:,iz,1] -
      (1.0 / fac) *AI.A23[:,:,iz,1] * AI.A32[:,:,iz,1]
    B = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-1)*6:iz*6]
    @views @. AI.B13[:,:,iz,1] = B[1:M2,5:6]
    @views @. AI.B23[:,:,iz,1] = B[M2+1:2*M2,5:6]
    @views @. AI.B32[:,:,iz,1] = B[2*M2+1:3*M2,3:4]
    C = JacP21[1+(iz-1)*6:iz*6,1+(iz-1)*M23:iz*M23]
    @views @. AI.C13[:,:,iz,1] = C[1:2,2*M2+1:3*M2]
    @views @. AI.C23[:,:,iz,1] = C[3:4,2*M2+1:3*M2]
    @views @. AI.C32[:,:,iz,1] = C[5:6,M2+1:2*M2]
  end
end


function BlockGauss(A11,A12,A21,A22,r1,r2)
  A11F = collect(A11)
  A12F = collect(A12)
  A21F = collect(A21)
  A22F = collect(A22)
  S = A22F - A21F * (A11F \ A12F)
  rS = r2 - (A21F * (A11F \ r1))
  x2 = S \ rS
  x1 = A11F \ (r1 - A12F * x2)
  SS = sparse(S)
  return x1, x2, SS
end

function SchurBoundary!(AI,A22)
  M2 = AI.M - 2
  invfac = 1.0 / AI.fac
  FacGrav = AI.FacGrav
  A22F = collect(A22)
  r1 = zeros(M2,2)
  r2 = zeros(M2,2)
  r3 = zeros(M2,2)
  s1 = zeros(2,2)
  s2 = zeros(2,2)
  s3 = zeros(2,2)
  for iz = 1 : nz
    sh = (iz - 1) * 6  
#   B32 3:4  
    @. r3 = -AI.B32[:,:,iz,1]   
    r3 = AI.SA[:,:,iz,1] \ r3
    r1 = -invfac * (AI.A13[:,:,iz,1] * r3)
    r2 = -invfac * (AI.A23[:,:,iz,1] * r3)
    s1 = -AI.C13[:,:,iz,1] * r3
    s2 = -AI.C23[:,:,iz,1] * r3
    @show size(AI.C32[:,:,iz,1])
    @show size(r2),size(r1)
    s3 = -AI.C32[:,:,iz,1] * r2  - FacGrav * sum(r1,dims=1)
    A22F[sh + 1:sh + 6,sh + 3:sh + 4] += [s1;s2;s3]
  end    
end

function Permutation(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  for iz = 1 : nz
    for iv = 1 : 3
      for k = 2 : M - 1
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
  end

#=
  iv = 1      
  for iz = 1 : nz
    ii += 1
    p[ii] = 1 + (iz-1) * M  + (iv - 1) * N
    ii += 1
    p[ii] = M + (iz - 1) * M + (iv - 1) * N
  end
=#  
  for iz = 1 : nz
    for iv = 1 : 3
      ii += 1
      p[ii] = 1 + (iz-1) * M  + (iv - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (iv - 1) * N
    end
  end
  return p
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

Problem = "HillAgnesiXCart"
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
Lx = 1000.0
Ly = 1000.0
x0 = 0.0
y0 = 0.0
nz = 20
OrdPoly = 4
OrdPolyZ = 6
M = OrdPolyZ + 1
N = M * nz
Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)

OrdPrint = 4
OrdPrintZ = 4
DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)
DG.NumG = 1
DG.NumI = 1

dSdS,dSdM,dMdS,dMdM = DGSEM.InitJacDG(DG,nz,Param)

fac = 0.5
U = ones(M,nz,DG.NumG,5)
H = 10000.0
dzLoc = H / nz
dz = ones(nz,DG.NumG) * dzLoc

Profile = Examples.StratifiedExample()(Param,Phys)
time = 0.0
for iz = 1 : nz
  z0 = (iz - 1) * dzLoc  
  z1 = iz * dzLoc  
  for k = 1 : M
    zLoc = 0.5 * ((1.0 - DG.xwZ[k]) * z0 + (1.0 + DG.xwZ[k]) * z1)  
    xS = SVector{3}(0.0,0.0,zLoc)
    RhoP,_,_,_,ThP= Profile(xS,time)
    U[k,iz,1,1] = RhoP
    U[k,iz,1,5] = RhoP * ThP
  end
end

r = ones(3*N)

JacLU, Jac  = DGSEM.JacDGT(U,DG,fac,dSdS,dSdM,dMdS,dMdM,dz,Phys)

p = Permutation(M,nz)

JacP = Jac[p,p]    

n1 = nz * (M-2) * 3

JacP11 = JacP[1:n1,1:n1]
JacP12 = JacP[1:n1,n1+1:end]
JacP21 = JacP[n1+1:end,1:n1]
JacP22 = JacP[n1+1:end,n1+1:end]

r = zeros(3*N)
aa = collect(1:3*N)
@. r = aa
rP = r[p]
rP1 = rP[1:n1]
rP2 = rP[n1+1:end]

x1P, x2P, SS = BlockGauss(JacP11,JacP12,JacP21,JacP22,rP1,rP2)
xP = [x1P; x2P]
xxS = xP[invperm(p)]
ldiv!(JacLU[1],r)
@show sum(abs.(xxS-r))

n1 = nz * 2
SS11 = SS[1:n1,1:n1]
SS12 = SS[1:n1,n1+1:end]
SS21 = SS[n1+1:end,1:n1]
SS22 = SS[n1+1:end,n1+1:end]
rP21 = rP2[1:n1]
rP22 = rP2[n1+1:end]

x21P, x22P, SSS = BlockGauss(SS11,SS12,SS21,SS22,rP21,rP22)

AI = Interior(M,nz,1)
FillInterior!(AI,JacP11,fac,Phys)
SchurBoundary!(AI,JacP22)
#=

#Individual inner block
M2 = M - 2
A13 = zeros(M2,M2,nz)
A23 = zeros(M2,M2,nz)
A31 = zeros(M2,M2,nz)
A32 = zeros(M2,M2,nz)
B13 = zeros(M2,2,nz)
B23 = zeros(M2,2,nz)
B32 = zeros(M2,2,nz)
C13 = zeros(2,M2,nz)
C23 = zeros(2,M2,nz)
C32 = zeros(2,M2,nz)
SA = zeros(M2,M2,nz)
for iz = 1 : nz
  A = JacP11[1+(iz-1)*M23:iz*M23,1+(iz-1)*M23:iz*M23]
  @views @. AI.A13[:,:,iz,1] = A[1:M2,2*M2+1:3*M2]
  @views @. AI.A23[:,:,iz,1] = A[M2+1:2*M2,2*M2+1:3*M2]
  @views @. AI.A31[:,:,iz,1] = A[2*M2+1:3*M2,1:M2]
  @views @. AI.A32[:,:,iz,1] = A[2*M2+1:3*M2,M2+1:2*M2]
  SA[:,:,iz] = fac * I - (0.5 * Phys.Grav / fac) * A13[:,:,iz] - (1.0 / fac) *A23[:,:,iz] * A32[:,:,iz]
  B = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-1)*6:iz*6]
  @views @. AI.B13[:,:,iz,1] = B[1:M2,5:6]
  @views @. AI.B23[:,:,iz,1] = B[M2+1:2*M2,5:6]
  @views @. AI.B32[:,:,iz,1] = B[2*M2+1:3*M2,3:4]
  C = JacP21[1+(iz-1)*6:iz*6,1+(iz-1)*M23:iz*M23]
  @views @. AI.C13[:,:,iz,1] = C[1:2,2*M2+1:3*M2]
  @views @. AI.C23[:,:,iz,1] = C[3:4,2*M2+1:3*M2]
  @views @. AI.C32[:,:,iz,1] = C[5:6,M2+1:2*M2]
end  
=#

iz = 2
M23 = (M-2) * 3
AA = JacP11[1+(iz-1)*M23:iz*M23,1+(iz-1)*M23:iz*M23]
CC = JacP21[1+(iz-1)*6:iz*6,1+(iz-1)*M23:iz*M23]
BB = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-1)*6:iz*6]

aa = 3

#=
A11 = A[1:M2,1:M2]
A12 = A[1:M2,M2+1:2*M2]
A13 = A[1:M2,2*M2+1:3*M2]
A21 = A[M2+1:2*M2,1:M2]
A22 = A[M2+1:2*M2,M2+1:2*M2]
A23 = A[M2+1:2*M2,2*M2+1:3*M2]
A31 = A[2*M2+1:3*M2,1:M2]
A32 = A[2*M2+1:3*M2,M2+1:2*M2]
A33 = A[2*M2+1:3*M2,2*M2+1:3*M2]
=#

