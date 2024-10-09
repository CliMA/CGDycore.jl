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

TestIter = 2

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

JuliaDevice = get(ENV, "JuliaDevice", "CPU")
JuliaGPU = get(ENV, "JuliaGPU", "CUDA")
machine = get(ENV, "machine", "")

if JuliaDevice == "CPU"
  backend = CPU()
elseif JuliaDevice == "GPU"
  if JuliaGPU == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(false)
    if machine == "levante"
    else
       CUDA.device!(Proc-1)
    end
  elseif JuliaGPU == "AMD"
    backend = ROCBackend()
    AMDGPU.allowscalar(false)
  elseif JuliaGPU == "Metal"
    backend = MetalBackend()
    Metal.allowscalar(true)
  end
else
  backend = CPU()
end




NzG = min(div(NumberThreadGPU,DoF),Nz)
group = (Ord, Ord, NzG, 1)
ndrange = (Ord, Ord, Nz, NF)
NumGG = min(div(NumberThreadGPU,Nz),NumG)
groupG = (Nz, NumGG)
ndrangeG = (Nz, NumG) 

F = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
@views FTh = F[:,:,5]
CacheF = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
U = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
@. U = abs(rand()) + 1
@views Th = U[:,:,5]
@views Rho = U[:,:,1]
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
KV = KernelAbstractions.ones(backend,FT,Nz,NumG)
KV = abs(rand()) + 1
dz = KernelAbstractions.ones(backend,FT,Nz,NumG)
dz = abs(rand()) + 100

KMomentumVectorInvariantCoriolisKernel! = GPU.MomentumVectorInvariantCoriolisKernel!(backend,group)
KMomentumVectorInvariantCoriolisKernel!(F,U,D,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Momentum"
@time for iter = 1 : TestIter
  KMomentumVectorInvariantCoriolisKernel!(F,U,D,dXdxI,J,X,M,Glob,CoriolisFun,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  

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
@show "UpwindRho and RhoTh"
@time for iter = 1 : TestIter
  KDivRhoThUpwind3Kernel!(F,U,D,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  

KDivRhoTrUpwind3Kernel! = GPU.DivRhoTrUpwind3Kernel!(backend,group)
KDivRhoTrUpwind3Kernel!(FTh,Th,U,D,dXdxI,J,M,Glob,ndrange=ndrange)
KernelAbstractions.synchronize(backend)
@show "Upwind Tracer"
@time for iter = 1 : TestIter
  KDivRhoTrUpwind3Kernel!(FTh,Th,U,D,dXdxI,J,M,Glob,ndrange=ndrange)
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

@show "Diffusion Scalar"
@. F = 0
KVerticalDiffusionScalarKernel! = GPU.VerticalDiffusionScalarKernel!(backend,groupG)
KVerticalDiffusionScalarKernel!(FTh,Th,Rho,KV,dz,ndrange=ndrangeG)
KernelAbstractions.synchronize(backend)
@show sum(abs.(F))
@time for iter = 1 : TestIter
  KVerticalDiffusionScalarKernel!(FTh,Th,Rho,KV,dz,ndrange=ndrangeG)
  KernelAbstractions.synchronize(backend)
end  

@show "DiffusionNew Scalar"
@. F = 0
KVerticalDiffusionScalarNewKernel! = GPU.VerticalDiffusionScalarNewKernel!(backend,groupG)
KVerticalDiffusionScalarNewKernel!(FTh,Th,Rho,KV,dz,ndrange=ndrangeG)
KernelAbstractions.synchronize(backend)
@show sum(abs.(F))
@time for iter = 1 : TestIter
  KVerticalDiffusionScalarNewKernel!(FTh,Th,Rho,KV,dz,ndrange=ndrangeG)
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
@show "MicroPhysics"
@time for iter = 1 : TestIter
  KMicrophysicsKernel!(MicrophysicsSource,F,U,p,ndrange=ndrangeG)
  KernelAbstractions.synchronize(backend)
end  


Param = Examples.Parameters(FT,"HeldSuarezMoistSphere")
_, Force = Examples.HeldSuarezDryExample()(Param,Phys)
KForceKernel! = GPU.ForceKernel!(backend, groupG)
KForceKernel!(Force,F,U,p,xS,ndrange=ndrangeG)
KernelAbstractions.synchronize(backend)
@show "HeldSuarez"
@time for iter = 1 : TestIter
  KForceKernel!(Force,F,U,p,xS,ndrange=ndrangeG)
  KernelAbstractions.synchronize(backend)
end  
@show "Ende"
