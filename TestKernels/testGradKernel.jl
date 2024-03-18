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

NF = 5400 # 30*30*6
NumG = 48602 
Ord = 4
DoF = Ord * Ord
NumV = 5
NumTr = 2
Nz = 64
NumberThreadGPU = 1024
GlobCPU = zeros(Int,DoF,NF)
read!("TestAMD/GlobInd",GlobCPU)



backend = CPU()
#backend = CUDABackend()
#backend = ROCBackend()
NzG = min(div(NumberThreadGPU,DoF),Nz)
group = (Ord, Ord, NzG, 1)
ndrange = (Ord, Ord, Nz, NF)

F = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
CacheF = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
U = KernelAbstractions.ones(backend,FT,Nz,NumG,NumV + NumTr)
D = KernelAbstractions.ones(backend,FT,4,4)
DW = KernelAbstractions.ones(backend,FT,4,4)
dXdxI = KernelAbstractions.ones(backend,FT,3,3,2,DoF,Nz,NF)
J = KernelAbstractions.ones(backend,FT,DoF,2,Nz,NF)
X = KernelAbstractions.ones(backend,FT,DoF,2,3,Nz,NF)
MRho = KernelAbstractions.ones(backend,FT,Nz-1,NumG)
M = KernelAbstractions.ones(backend,FT,Nz,NumG)
p = KernelAbstractions.ones(backend,FT,Nz,NumG)
KoeffCurl = FT(1)
KoeffGrad = FT(1)
KoeffDiv = FT(1)
Glob = KernelAbstractions.zeros(backend,Int,DoF,NF)
copyto!(Glob,GlobCPU)

KGradKernel! = GPU.GradKernel!(backend,group)
KGradKernel!(F,U,p,D,dXdxI,J,M,MRho,Glob,Phys,ndrange=ndrange)
KernelAbstractions.synchronize(backend)

