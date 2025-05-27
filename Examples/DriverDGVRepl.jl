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


nz = 4
H = 10000.0
OrdPolyZ = 4
OrdPrintZ = 4

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
UStart = similar(U)
UImp = similar(U)
UStart .= U

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
FluxA = DGVertical.EulerFluxA()(RhoPos,wPos,ThPos,pAuxPos)
Flux = DGVertical.EulerFlux()(RhoPos,wPos,ThPos,pAuxPos)
RiemannSolver = DGVertical.RiemannLMARSV()(Param,Phys,RhoPos,wPos,ThPos,pAuxPos)
RiemannSolverA = DGVertical.RiemannLMARSAV()(Param,Phys,RhoPos,wPos,ThPos,pAuxPos)

N = nz*(OrdPolyZ+1)

dSdS,dSdM = DGVertical.DScalarDMomAc(nz,DG1,Param.cS)
dMdS,dMdM = DGVertical.DMomDScalarAc(nz,DG1,Param.cS)

dSdS = spdiagm(reshape(1.0./J,N)) * dSdS
dSdM = spdiagm(reshape(1.0./J,N)) * dSdM
dMdS = spdiagm(reshape(1.0./J,N)) * dMdS
dMdM = spdiagm(reshape(1.0./J,N)) * dMdM

# S => Th
# M => w

#=
Th = reshape(U[:,:,ThPos]./U[:,:,RhoPos],N)
dpdRhoTh = reshape(Phys.Rd * (Phys.Rd * U[:,:,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
spId = sparse(I,N,N)
spZero = spzeros(N,N)

JacAc = [spZero dMdS * diagm(reshape(dpdRhoTh,N))
         dSdM* diagm(reshape(Th,N))  spZero]
JacAcFull = zeros(2*N,2*N)
JacAcFull .= JacAc

Jac = [spzeros(N,N) dSdM spzeros(N,N)
       -Phys.Grav*spId spZero dMdS * diagm(reshape(dpdRhoTh,N))
       spZero dSdM* diagm(reshape(Th,N))  spZero]
JacFull = zeros(3*N,3*N)
JacFull .= Jac

stop         
=#
Th = reshape(U[:,:,ThPos]./U[:,:,RhoPos],N)
@. Th = 300.0
dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd * 
  (Phys.Rd * U[:,:,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
JacStart = [spzeros(N,N) dSdM  dSdS* diagm(dpdRhoTh)
       -Phys.Grav * sparse(I,N,N) dMdM  dMdS * diagm(dpdRhoTh)
       spzeros(N,N) dSdM* diagm(Th)  diagm(Th) * dSdS* diagm(dpdRhoTh)]
JacDiff = zeros(3*N,3*N)
ColJacDiff = zeros(OrdPolyZ + 1,nz,3)
DGVertical.FcnGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxA,RiemannSolverA) 
iC = 0
for iv = 1 : 3
  for iz = 1 : nz
    for k = 1 : OrdPolyZ + 1
      global iC += 1  
      temp = U[k,iz,iv]
      U[k,iz,iv] = (1 + 1.e-8) * U[k,iz,iv] + 1.e-8
      DGVertical.FcnGPUVert!(FT,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxA,RiemannSolverA) 
      @. ColJacDiff .= (FT - F) / (U[k,iz,iv] - temp)
      JacDiff[:,iC] .= reshape(ColJacDiff,3*N)
      U[k,iz,iv] = temp
    end
  end
end  
stop

    

dtau = 50.0
fac = dtau 

nIter = 100
@. U = UStart
for Iter = 1 : nIter
   @show "N",Iter,sum(abs.(U[:,:,wPos]))   
   DGVertical.FcnGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,Flux,RiemannSolver) 
   Th = reshape(U[:,:,ThPos]./U[:,:,RhoPos],N)
   dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd * 
     (Phys.Rd * U[:,:,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
   Jac = [spzeros(N,N) dSdM  dSdS* diagm(dpdRhoTh)
       -Phys.Grav * sparse(I,N,N) dMdM  dMdS * diagm(dpdRhoTh)
       spzeros(N,N) dSdM* diagm(Th)  diagm(Th) * dSdS* diagm(dpdRhoTh)]   
   A = (1 / fac) * sparse(I,3*N,3*N) - Jac       
   k = reshape(A \ reshape(F,3*N),OrdPolyZ+1,nz,3)
   @. U  = U + k
end   
UImp .= U

@. U = UStart
for Iter = 1 : nIter
   @show "S",Iter,sum(abs.(U[:,:,wPos]))
   DGVertical.FcnSplitGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver)
   Th = reshape(U[:,:,ThPos]./U[:,:,RhoPos],N)
   dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
     (Phys.Rd * U[:,:,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
   Jac = [spzeros(N,N) dSdM  dSdS* diagm(dpdRhoTh)
       -Phys.Grav * sparse(I,N,N) dMdM  dMdS * diagm(dpdRhoTh)
       spzeros(N,N) dSdM* diagm(Th)  diagm(Th) * dSdS* diagm(dpdRhoTh)]
   A = (1 / fac) * sparse(I,3*N,3*N) - Jac
   k = reshape(A \ reshape(F,3*N),OrdPolyZ+1,nz,3)
   @. U  = U + k
end
UImp .= U
stop

dtau = .2
nIter = 1000
U .= UStart
for Iter = 1 : nIter
   DGVertical.FcnGPUVert!(F,U,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/3 * dtau *F
   DGVertical.FcnGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. UNew = U + 1/2 * dtau *F
   DGVertical.FcnGPUVert!(F,UNew,DG1,X,dXdxI,J,CacheU,Pressure,Phys,FluxAverage,RiemannSolver) 
   @. U = U + dtau *F
 end  
 @show "Ende Exp",sum(abs.(U[:,:,wPos]))   
