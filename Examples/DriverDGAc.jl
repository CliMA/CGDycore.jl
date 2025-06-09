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


nz = 10
H = 10000.0
OrdPolyZ = 5
OrdPrintZ = 4

z,zP,dzeta = Grids.AddVerticalGrid(nz,H)

DG1 = FiniteElements.DG1{FTB}(backend,OrdPolyZ,OrdPrintZ)

X = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)
dXdxI = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)
J = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz)

Grids.JacobiDG1GPU!(X,dXdxI,J,DG1,z)
NumV = 2
NumAux = 0
UStart = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz,NumV)
Profile = Examples.StratifiedExample()(Param,Phys)
time = 0.0
for iz = 1 : nz
  for K = 1 : DG1.OrdPolyZ + 1
    xS = SVector{3}(0.0,0.0,X[K,iz])
    UStart[K,iz,1] = cos(pi * X[K,iz] / H)
    UStart[K,iz,2] = sin(pi * X[K,iz] / H)
  end
end  

F = similar(UStart)
FT = similar(UStart)
U = similar(UStart)
UNew = similar(UStart)
ULin = similar(UStart)
UImp = similar(UStart)
USplit = similar(UStart)

CacheU = KernelAbstractions.zeros(backend,FTB,DG1.OrdPolyZ+1,nz,NumAux)
Pressure, dPresdRhoTh, dPresdRho = Models.DryDG()(Phys)

pPos = 1
wPos = 2

cS = 360.0
FluxAverage = DGVertical.KennedyGruberAccousticV()(pPos,wPos,cS)
Flux = DGVertical.AccousticFlux()(pPos,wPos,cS)
RiemannSolver = DGVertical.RiemannAccousticV()(Param,Phys,pPos,wPos,cS)

N = nz*(OrdPolyZ+1)
dSdS,dSdM = DGVertical.DScalarDMomAc(nz,DG1,cS)
dSdS = spdiagm(reshape(1.0./J,N)) * dSdS
dSdM = spdiagm(reshape(1.0./J,N)) * dSdM

dMdS,dMdM = DGVertical.DMomDScalarAc(nz,DG1,cS)
dMdS = spdiagm(reshape(1.0./J,N)) * dMdS
dMdM = spdiagm(reshape(1.0./J,N)) * dMdM

Jac = [cS^2 * dSdS cS^2 * dSdM
     dMdS dMdM]
DGVertical.FcnAccousticGPUVert!(F,UStart,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
JacDiff = zeros(2*N,2*N)
ColJacDiff = zeros(OrdPolyZ+1,nz,2)
iC = 0
for iv = 1 : 2
  for iz = 1 : nz
    for k = 1 : OrdPolyZ + 1
      global iC += 1
      temp = UStart[k,iz,iv]
      UStart[k,iz,iv] = (1 + 1.e-8) * UStart[k,iz,iv] + 1.e-8
      DGVertical.FcnAccousticGPUVert!(FT,UStart,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
      ColJacDiff .= (FT - F) ./ (UStart[k,iz,iv] - temp)
      JacDiff[:,iC] .= reshape(ColJacDiff,2*N)
      UStart[k,iz,iv] = temp
    end
  end
end
aa = Jac * reshape(UStart,2*N)
aa = reshape(aa,OrdPolyZ+1,nz,2)
@show sum(abs.(aa-F))

nIter = 1000
dtau = 0.1 * H / N / cS
@. U = UStart
for Iter = 1 : nIter
#  DGVertical.FcnAccousticGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
   F = reshape(Jac * reshape(U,2*N),OrdPolyZ+1,nz,2)
   @. UNew = U + 1/3 * dtau *F
#  DGVertical.FcnAccousticGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
   F = reshape(Jac * reshape(UNew,2*N),OrdPolyZ+1,nz,2)
   @. UNew = U + 1/2 * dtau *F
#  DGVertical.FcnAccousticGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
   F = reshape(Jac * reshape(UNew,2*N),OrdPolyZ+1,nz,2)
   @. U = U + dtau *F
end
@. ULin = U



nIter = 100
dtau = H / N / cS
fac = dtau 
Jac = 1/fac * sparse(I,2*N,2*N) - Jac
@. U = UStart
for Iter = 1 : nIter
   DGVertical.FcnSplitAccousticGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   k = reshape(Jac \ reshape(F,NumV*N),OrdPolyZ+1,nz,NumV)
   @. U  = U + k
   @show Iter,U[:,1,wPos]
end   
@. UImp = U

nIter = 1000
dtau = 0.1 * H / N / cS
@. U = UStart
for Iter = 1 : nIter
   @show sum(abs.(U))  
   DGVertical.FcnSplitAccousticGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/3 * dtau *F
   DGVertical.FcnSplitAccousticGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/2 * dtau *F
   DGVertical.FcnSplitAccousticGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. U = U + dtau *F
end  
@. USplit = U

nIter = 1000
dtau = 0.1 * H / N / cS
@. U = UStart
for Iter = 1 : nIter
   DGVertical.FcnAccousticGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
   @. UNew = U + 1/3 * dtau *F
   DGVertical.FcnAccousticGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
   @. UNew = U + 1/2 * dtau *F
   DGVertical.FcnAccousticGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver)
   @. U = U + dtau *F
end
