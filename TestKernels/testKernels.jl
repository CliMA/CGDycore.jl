import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using MPI


FT = Float32
Phys = DyCore.PhysParameters{FT}()

TestIter = 200

NF = 5400 # 30*30*6
NumG = 48602 
Ord = 4
DoF = Ord * Ord
NumV = 5
NumTr = 2
Nz = 64
NumberThreadGPU = 512
GlobCPU = zeros(Int,DoF,NF)

read!("TestKernels/GlobInd",GlobCPU)



#backend = CPU()
backend = CUDABackend()
#backend = ROCBackend()
NzG = min(div(NumberThreadGPU,DoF),Nz)
group = (Ord, Ord, NzG, 1)
ndrange = (Ord, Ord, Nz, NF)
NumGG = min(div(NumberThreadGPU,Nz),NumG)
groupG = (Nz, NumGG)
ndrangeG = (Nz, NumG) 

F = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
CacheF = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
U = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
D = KernelAbstractions.ones(backend,FT,4,4)
DW = KernelAbstractions.ones(backend,FT,4,4)
dXdxI = KernelAbstractions.ones(backend,FT,3,3,2,DoF,Nz,NF)
J = KernelAbstractions.ones(backend,FT,DoF,2,Nz,NF)
X = KernelAbstractions.ones(backend,FT,DoF,2,3,Nz,NF)
xS    = KernelAbstractions.zeros(backend,FT,2,NumG)
M = KernelAbstractions.ones(backend,FT,Nz,NumG,2)
p = KernelAbstractions.ones(backend,FT,Nz,NumG)
KoeffCurl = FT(1)
KoeffGrad = FT(1)
KoeffDiv = FT(1)
Glob = KernelAbstractions.zeros(backend,Int,DoF,NF)
copyto!(Glob,GlobCPU)
CoriolisFun = GPU.CoriolisShallow()(Phys)
GravitationFun = GPU.GravitationShallow()(Phys)

KMomentumVectorInvariantCoriolisKernel! = GPU.MomentumVectorInvariantCoriolisKernel!(backend,group)
KMomentumVectorInvariantCoriolisKernel!(F,U,D,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Momentum"
@time for iter = 1 : TestIter
  KMomentumVectorInvariantCoriolisKernel!(F,U,D,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  
@show CUDA.@profile trace=true CUDA.@sync KMomentumVectorInvariantCoriolisKernel!(F,U,D,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrange)
KernelAbstractions.synchronize(backend)

KGradKernel! = GPU.GradKernel!(backend,group)
KGradKernel!(F,U,p,D,dXdxI,J,X,M,Glob,GravitationFun,ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Gradient"
@time for iter = 1 : TestIter
  KGradKernel!(F,U,p,D,dXdxI,J,X,M,Glob,GravitationFun,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  

KDivRhoThUpwind3Kernel! = GPU.DivRhoThUpwind3Kernel!(backend,group)
KDivRhoThUpwind3Kernel!(F,U,D,dXdxI,J,M,Glob,ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Upwind"
@time for iter = 1 : TestIter
  KDivRhoThUpwind3Kernel!(F,U,D,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  

KHyperViscKoeffKernel! = GPU.HyperViscKoeffKernel!(backend,group)
KHyperViscKoeffKernel!(F,U,CacheF,D,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "HyperVisc"
@time for iter = 1 : TestIter
  KHyperViscKoeffKernel!(F,U,CacheF,D,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  

RhoPos = 1
ThPos = 5
RhoVPos = 6
RhoCPos = 7
RelCloud = 1.e-1
Rain = 0.0
MicrophysicsSource  = Models.SimpleMicrophysics()(Phys,RhoPos,ThPos,
RhoVPos,RhoCPos,RelCloud,Rain)
KMicrophysicsKernel! = GPU.MicrophysicsKernel!(backend, groupG)
KMicrophysicsKernel!(MicrophysicsSource,F,U,p,ndrange=ndrangeG)
KernelAbstractions.synchronize(backend)
@time for iter = 1 : TestIter
  KMicrophysicsKernel!(MicrophysicsSource,F,U,p,ndrange=ndrangeG)
  KernelAbstractions.synchronize(backend)
end  


Param = Examples.Parameters(FT,"HeldSuarezMoistSphere")
_, Force = Examples.HeldSuarezDryExample()(Param,Phys)
KForceKernel! = GPU.ForceKernel!(backend, groupG)
KForceKernel!(Force,F,U,p,xS,ndrange=ndrangeG)
KernelAbstractions.synchronize(backend)
@time for iter = 1 : TestIter
  KForceKernel!(Force,F,U,p,xS,ndrange=ndrangeG)
  KernelAbstractions.synchronize(backend)
end  
stop
