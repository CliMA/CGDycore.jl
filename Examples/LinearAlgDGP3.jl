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
using BandedMatrices, CliqueTrees

mutable struct Interior 
  M::Int
  nz::Int
  fac::Float64
  FacGrav::Float64
  A13::Array{Float64, 4}
  A23::Array{Float64, 4}
  A32::Array{Float64, 4}
  BOD::Array{Float64, 3}
  B11::Array{Float64, 2}
  B12::Array{Float64, 4}
  B21::Array{Float64, 2}
  B22::Array{Float64, 4}
  B31::Array{Float64, 4}
  C13::Array{Float64, 4}
  C22::Array{Float64, 4}
  luSA::Array{LU, 2}
end  

function Interior(M,nz,NumG)
  M2 = M - 2
  fac = 0
  FacGrav = 0
  A13 = zeros(M,M2,nz,NumG)
  A23 = zeros(M2,M2,nz,NumG)
  A32 = zeros(M2,M2,nz,NumG)
  BOD = zeros(4,nz,NumG)
  B11 = zeros(nz,NumG)
  B12 = zeros(M,2,nz,NumG)
  B21 = zeros(nz,NumG)
  B22 = zeros(M2,2,nz,NumG)
  B31 = zeros(M2,2,nz,NumG)
  C13 = zeros(2,M2,nz,NumG)
  C22 = zeros(2,M2,nz,NumG)
  luSA = Array{LU}(undef,nz,NumG)

  return Interior(
    M,
    nz,
    fac,
    FacGrav,
    A13,
    A23,
    A32,
    BOD,
    B11,
    B12,
    B21,
    B22,
    B31,
    C13,
    C22,
    luSA,
  )
end  

function FillInterior!(AI,JacP11,JacP12,JacP21,fac,Phys)
  M2 = AI.M - 2
  M23= M2 * 2 + M
  nz = AI.nz
  AI.fac = fac
  AI.FacGrav = 0.5 * Phys.Grav
  SA = zeros(M2,M2)

  for iz = 1 : nz
    A = JacP11[1+(iz-1)*M23:iz*M23,1+(iz-1)*M23:iz*M23]
    @views @. AI.A13[:,:,iz,1] = A[1:M,M+M2+1:end]
    @views @. AI.A23[:,:,iz,1] = A[M+1:M+M2,M+M2+1:end]
    @views @. AI.A32[:,:,iz,1] = A[M+M2+1:end,M+1:M+M2]
    @views SA = fac * I - (0.5 * Phys.Grav / fac) * AI.A13[2:M-1,:,iz,1] - 
      (1.0 / fac) *AI.A32[:,:,iz,1] * AI.A23[:,:,iz,1]
    AI.luSA[iz,1] = lu(SA)
    B = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-1)*4:iz*4]
    AI.B11[iz,1] = B[1,1]
    @views @. AI.B12[:,:,iz,1] = B[1:M,3:4]
    AI.B21[iz,1] = B[M,2]
    if iz > 1
      BM = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-2)*4:(iz-1)*4]  
      AI.BOD[1,iz,1] = BM[1,2]
      AI.BOD[2,iz,1] = BM[1,4]
    end  
    if iz < nz
      BP = JacP12[1+(iz-1)*M23:iz*M23,1+iz*4:(iz+1)*4]
      AI.BOD[3,iz,1] = BP[M,1]
      AI.BOD[4,iz,1] = BP[M,3]
    end  
    @views @. AI.B22[:,:,iz,1] = B[1+M:M+M2,3:4]
    @views @. AI.B31[:,:,iz,1] = B[M+M2+1:end,1:2]
    C = JacP21[1+(iz-1)*4:iz*4,1+(iz-1)*M23:iz*M23]
    @views @. AI.C13[:,:,iz,1] = C[1:2,M+M2+1:end]
#   C21 = FacGrav Position 1 and N
    @views @. AI.C13[:,:,iz,1] = C[1:2,M+M2+1:end]
    @views @. AI.C22[:,:,iz,1] = C[3:4,M+1:M+M2]
  end
end
function FillInteriorDirect!(AI,AID,JacP11,JacP12,JacP21,U,DG,dz,fac,Phys,Param)
  M = DG.OrdPolyZ + 1
  M1 = M - 1
  M2 = M - 2
  M23= M2 * 2 + M
  nz = size(dz,1)
  AI.fac = fac
  AI.FacGrav = 0.5 * Phys.Grav
  SA = zeros(M2,M2)
  DW = DG.DWZ
  ThPos = 5
  RhoPos = 1
  Th = zeros(M)
  dpdRhoTh = zeros(M)
  cS = Param.cS

  ID = 1
  for iz = 1 : nz
    Th = U[:,iz,ID,ThPos]./U[:,iz,ID,RhoPos]
    dpdRhoTh = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,iz,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
    A = JacP11[1+(iz-1)*M23:iz*M23,1+(iz-1)*M23:iz*M23]
    @views @. AI.A13[:,:,iz,1] = A[1:M,M+M2+1:end]
    @views @. AID.A13[:,:,iz,1] = (2 / dz[iz]) * DW[:,2:M1]

    @views @. AI.A23[:,:,iz,1] = A[M+1:M+M2,M+M2+1:end]
    @views AID.A23[:,:,iz,1] = (2 / dz[iz]) * DW[2:M1,2:M1] * diagm(Th[2:M-1])

    @views @. AI.A32[:,:,iz,1] = A[M+M2+1:end,M+1:M+M2]
    @views AID.A32[:,:,iz,1] = (2 / dz[iz]) * DW[2:M1,2:M1] * diagm(dpdRhoTh[2:M-1])

    @views SA = fac * I - (0.5 * Phys.Grav / fac) * AI.A13[2:M-1,:,iz,1] - 
      (1.0 / fac) *AI.A32[:,:,iz,1] * AI.A23[:,:,iz,1]
    @views SAD = fac * I - (0.5 * Phys.Grav / fac) * AID.A13[2:M-1,:,iz,1] - 
      (1.0 / fac) *AID.A32[:,:,iz,1] * AID.A23[:,:,iz,1]
    AI.luSA[iz,1] = lu(SA)
    AID.luSA[iz,1] = lu(SAD)

    B = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-1)*4:iz*4]

    AI.B11[iz,1] = B[1,1]
    AID.B11[iz,1] = -(1 / dz[iz]) * fac/cS/DG.wZ[1] * dpdRhoTh[1]
    @views @. AI.B12[:,:,iz,1] = B[1:M,3:4]
    AI.B21[iz,1] = B[M,2]
    if iz > 1
      BM = JacP12[1+(iz-1)*M23:iz*M23,1+(iz-2)*4:(iz-1)*4]  
      AI.BOD[1,iz,1] = BM[1,2]
      AI.BOD[2,iz,1] = BM[1,4]
    end  
    if iz < nz
      BP = JacP12[1+(iz-1)*M23:iz*M23,1+iz*4:(iz+1)*4]
      AI.BOD[3,iz,1] = BP[M,1]
      AI.BOD[4,iz,1] = BP[M,3]
    end  
    @views @. AI.B22[:,:,iz,1] = B[1+M:M+M2,3:4]
    @views @. AI.B31[:,:,iz,1] = B[M+M2+1:end,1:2]
    C = JacP21[1+(iz-1)*4:iz*4,1+(iz-1)*M23:iz*M23]
    @views @. AI.C13[:,:,iz,1] = C[1:2,M+M2+1:end]
#   C21 = FacGrav Position 1 and N
    @views @. AI.C13[:,:,iz,1] = C[1:2,M+M2+1:end]
    @views @. AI.C22[:,:,iz,1] = C[3:4,M+1:M+M2]
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

function SchurBoundary!(AI,A22B)
  M2 = AI.M - 2
  invfac = 1.0 / AI.fac
  FacGrav = AI.FacGrav
  r1 = zeros(M,2)
  r2 = zeros(M2,2)
  r3 = zeros(M2,2)
  s1 = zeros(2,2)
  s2 = zeros(2,2)
  rr1 = zeros(M,1)
  rr2 = zeros(M2,1)
  rr3 = zeros(M2,1)
  ss1 = zeros(2,1)
  ss2 = zeros(2,1)
  r11 = zeros(1,2)
  r11 = zeros(1,1)
  r1M = zeros(1,1)
  for iz = 1 : nz
    sh = (iz - 1) * 4  
#   Column 1:2  
    @. r3 = AI.B31[:,:,iz,1]   

    ldiv!(AI.luSA[iz,1],r3)

    r11 = -invfac * (AI.A13[1:1,:,iz,1] * r3[:,:])
    r1M = -invfac * (AI.A13[M:M,:,iz,1] * r3[:,:])
    r11[1,1] += invfac * AI.B11[iz,1]
    r1M[1,2] += invfac * AI.B21[iz,1]
    @views r2 = -invfac * (AI.A23[:,:,iz,1] * r3)

    @views s1 = -AI.C13[:,:,iz,1] * r3
    @views s2 = -AI.C22[:,:,iz,1] * r2
    s2[1,1] = s2[1,1] - FacGrav * r11[1,1]
    s2[2,1] = s2[2,1] - FacGrav * r1M[1,1]
    s2[1,2] = s2[1,2] - FacGrav * r11[1,2]
    s2[2,2] = s2[2,2] - FacGrav * r1M[1,2]
    A22B[sh + 1:sh + 4,sh + 1:sh + 2] .+= [s1;s2]

#   Column 3:4  
    @. r1 = AI.B12[:,:,iz,1]
    @. r2 = AI.B22[:,:,iz,1]
    @views r3 = -invfac * (AI.A32[:,:,iz,1] * r2 + FacGrav * r1[2:M-1,:]) 

    ldiv!(AI.luSA[iz,1],r3)

    r11 = invfac * (r1[1:1,:] - AI.A13[1:1,:,iz,1] * r3[:,:])
    r1M = invfac * (r1[M:M,:] - AI.A13[M:M,:,iz,1] * r3[:,:])
    @views r2 = invfac * (r2 - AI.A23[:,:,iz,1] * r3)

    @views s1 = -AI.C13[:,:,iz,1] * r3
    @views s2 = -AI.C22[:,:,iz,1] * r2
    s2[1,1] = s2[1,1] - FacGrav * r11[1,1]
    s2[2,1] = s2[2,1] - FacGrav * r1M[1,1]
    s2[1,2] = s2[1,2] - FacGrav * r11[1,2]
    s2[2,2] = s2[2,2] - FacGrav * r1M[1,2]
    A22B[sh + 1:sh + 4,sh + 3:sh + 4] .+= [s1;s2]

    @show iz,A22B[11,9]
    if iz > 1
#     Column -3  
#     @show "C-3",FacGrav * invfac * AI.BOD[1,iz,1]
      A22B[sh + 3,sh -2] -= FacGrav * invfac * AI.BOD[1,iz,1]

      #     Column -1
#     @show "C-1",FacGrav * invfac * AI.BOD[2,iz,1]
      A22B[sh + 3,sh] -= FacGrav * invfac * AI.BOD[2,iz,1]
    end
    if iz < nz 
#     Column +1  
#     @show "C+1",FacGrav * invfac * AI.BOD[3,iz,1]
      A22B[sh + 4,sh + 5] -= FacGrav * invfac * AI.BOD[3,iz,1]
#     Column +3  
#     @show "C+3",FacGrav * invfac * AI.BOD[4,iz,1]
      A22B[sh + 4,sh + 7] -= FacGrav * invfac * AI.BOD[4,iz,1]
    end
  end    
end

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
    for iv = 2 : 3
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
OrdPolyZ = 5
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
fac = 1.0
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

n1 = nz * (M-2) * 2 + nz * M

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
AID = Interior(M,nz,1)
FillInterior!(AI,JacP11,JacP12,JacP21,fac,Phys)
FillInteriorDirect!(AI,AID,JacP11,JacP12,JacP21,U,DG,dz,fac,Phys,Param)
JacP22B = BandedMatrix(JacP22)
SSB = BandedMatrix(SS)
luSSB = lu(SSB)
SchurBoundary!(AI,JacP22B)
stop
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

