import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGVertical, GPU, DyCore
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


backend = CPU()
FTB = Float64
Parallel = true

Problem = "HillAgnesiXCart"
Param = Examples.Parameters(FTB,Problem)
# Physical parameters
Phys = DyCore.PhysParameters{FTB}()


nz = 20
H = 10000.0
OrdPolyZ = 4
OrdPrintZ = 4
M = OrdPolyZ + 1
N = M * nz

z,zP,dzeta = Grids.AddVerticalGrid(nz,H)

DG1 = FiniteElements.DG1{FTB}(backend,OrdPolyZ,OrdPrintZ)

X = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)
dXdxI = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)
J = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)

Grids.JacobiDG1GPU!(X,dXdxI,J,DG1,z)
NumV = 3
NumAux = 2
U = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz,NumV)
Profile = Examples.StratifiedExample()(Param,Phys)
time = 0.0
for iz = 1 : nz
  for K = 1 : DG1.OrdPolyZ + 1
    xS = SVector{3}(0.0,0.0,X[K,iz])
    RhoP,_,_,_,ThP= Profile(xS,time)
    U[K,iz,1] = RhoP
    U[K,iz,3] = RhoP * ThP
  end
end  
UE = similar(U)
UE .= U
for iz = 1 : nz
  for K = 1 : DG1.OrdPolyZ + 1
    xS = SVector{3}(0.0,0.0,X[K,iz])
    if X[K,iz] > 500.0 && X[K,iz] <  2500.0
      U[K,iz,3] = U[K,iz,3] + 2.0
    end
  end
end
UStart = similar(U)
UImp1 = similar(U)
UImp2 = similar(U)
UStart .= U

@views @. UStart[:,3,3] += 3.0

F = similar(U)
FT = similar(U)
UNew = similar(U)
CacheU = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz,NumAux)
Pressure, dPresdRhoTh, dPresdRho = Models.DryDG()(Phys)

RhoPos = 1
wPos = 2
ThPos = 3
pAuxPos = 1
GPAuxPos = 2

FluxAverage = DGVertical.KennedyGruberGravV()(RhoPos,wPos,ThPos,pAuxPos,GPAuxPos)
Flux = DGVertical.EulerFlux()(RhoPos,wPos,ThPos,pAuxPos)
RiemannSolver = DGVertical.RiemannLMARSV()(Param,Phys,RhoPos,wPos,ThPos,pAuxPos)


dSdS,dSdM,dMdS,dMdM = DGVertical.InitJacDG(DG1,nz,Param)

#=
# Play ground for own linear algebra
Jac = DGVertical.JacDG1(U,DG1,dSdS,dSdM,dMdS,dMdM,J,Phys)
fac = 50.0
A = (1 / fac) * sparse(I,3*N,3*N) .- Jac[1]
 A = [A11 A12 A13
      A21 A22 A23
      A31 A32 A33]
 A = [A11 A13 A12
      A31 A33 A32
      A21 A23 A22]
 S = A22 - [A21 A23]*[A11 A13
                      A31 A33]^-1*[A12
                                   A32]
A11 = A[1:N,1:N]
A12 = A[1:N,N+1:2*N]
A13 = A[1:N,2*N+1:3*N]
A21 = A[N+1:2*N,1:N]
A22 = A[N+1:2*N,N+1:2*N]
A23 = A[N+1:2*N,2*N+1:3*N]
A31 = A[2*N+1:3*N,1:N]
A32 = A[2*N+1:3*N,N+1:2*N]
A33 = A[2*N+1:3*N,2*N+1:3*N]
D11 = [A11 A13
       A31 A33]
D12 = [A12
       A32]
D21 = [A21 A23]       
D22 = A22

DT = [A22 A23
      A32 A33]

D12F = zeros(2*N,N)
D12F .= D12
SchurF = D22 - D21 * (D11 \ D12F)
Schur = sparse(SchurF)
permuS = zeros(Int,N)
iInd = 1
for iz = 1 : nz
  for i = 2 : M - 1
    permuS[iInd] = i + (iz-1)*M
    global iInd += 1
  end
end
for iz = 1 : nz
  permuS[iInd] = 1 + (iz-1)*M
  global iInd += 1
  permuS[iInd] = M + (iz-1)*M
  global iInd += 1
end
SchurP = Schur[permuS,permuS]


permu = zeros(Int,3*N)
iInd = 1
  for iz = 1 : nz
    for i = 2 : M - 1
for iv = 1 : NumV
      permu[iInd] = i + (iz-1)*M + (iv-1)*M*nz
      global iInd += 1
    end
  end
end
  for iz = 1 : nz
for iv = 1 : NumV
    permu[iInd] = 1 + (iz-1)*M + (iv-1)*M*nz
    global iInd += 1
    permu[iInd] = M + (iz-1)*M + (iv-1)*M*nz
    global iInd += 1
  end
end
AP = A[permu,permu]
=#

dtau = 10.0
fac = dtau 

nIter = 100
@. U = UStart
for Iter = 1 : nIter
   @show "N",Iter,sum(abs.(U[:,:,wPos]))   
   DGVertical.FcnGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver) 
   Jac = DGVertical.JacDG1(U,DG1,dSdS,dSdM,dMdS,dMdM,J,Phys)
   A = (1 / fac) * sparse(I,3*N,3*N) - Jac[1]       
   k = reshape(A \ reshape(F,3*N),OrdPolyZ+1,nz,3)
   @. U  = U + k
end   
UImp1 .= U

@. U = UStart
for Iter = 1 : nIter
   @show "S",Iter,sum(abs.(U[:,:,wPos]))
   DGVertical.FcnSplitGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver)
   Jac = DGVertical.JacDG1(U,DG1,dSdS,dSdM,dMdS,dMdM,J,Phys)
   A = (1 / fac) * sparse(I,3*N,3*N) - Jac[1]
   k = reshape(A \ reshape(F,3*N),OrdPolyZ+1,nz,3)
   @. U  = U + k
end
UImp2 .= U

dtau = .1
nIter = 1000
U .= UStart
for Iter = 1 : nIter
   DGVertical.FcnSplitGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/3 * dtau *F
   DGVertical.FcnSplitGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/2 * dtau *F
   DGVertical.FcnSplitGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. U = U + dtau *F
end  
@show "Ende Exp",sum(abs.(U[:,:,wPos]))   
