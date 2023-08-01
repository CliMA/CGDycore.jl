using CGDycore
using MPI
using Random
using CUDA
using CUDA.CUDAKernels
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace


function DivRhoTrColumnE!(FRhoTrC,uC,vC,w,RhoTrC,D,dXdxI)

  N = size(uC,1)
  Nz = size(uC,3)

  @inbounds for iz = 1 : Nz
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N
        temp1 = -RhoTrC[i,j,iz] * (dXdxI[i,j,1,iz,3,1] * uC[i,j,iz] +
          dXdxI[i,j,1,iz,3,2] * vC[i,j,iz] + dXdxI[i,j,1,iz,3,3] * w[i,j,iz])
        temp2 = -RhoTrC[i,j,iz] * (dXdxI[i,j,2,iz,3,1] * uC[i,j,iz] +
          dXdxI[i,j,2,iz,3,2] * vC[i,j,iz] + dXdxI[i,j,2,iz,3,3] * w[i,j,iz+1])
        FRhoTrC[i,j,iz] += (temp2 - temp1)
        if iz > 1
          FRhoTrC[i,j,iz-1] -= 1/2 * temp1
          FRhoTrC[i,j,iz] += 1/2 * temp2
        end
        if iz < Nz
          FRhoTrC[i,j,iz] -= 1/2 * temp2
          FRhoTrC[i,j,iz+1] += 1/2 * temp1
        end
      end
    end
  end
  @inbounds for iz = 1 : Nz
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N
        tempx = -RhoTrC[i,j,iz] * ((dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * uC[i,j,iz] +
          (dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * vC[i,j,iz] +
          dXdxI[i,j,1,iz,1,3] * w[i,j,iz] + dXdxI[i,j,2,iz,1,3] * w[i,j,iz+1])
        tempy = -RhoTrC[i,j,iz] * ((dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * uC[i,j,iz] +
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * vC[i,j,iz] +
          dXdxI[i,j,1,iz,2,3] * w[i,j,iz] + dXdxI[i,j,2,iz,2,3] * w[i,j,iz+1])
        for k = 1 : N
          FRhoTrC[k,j,iz] += D[k,i] * tempx
          FRhoTrC[i,k,iz] += D[k,j] * tempy
        end
      end
    end
  end
end

function DivRhoTrColumnE2!(FRhoTrC,uC,vC,w,RhoTrC,D,dXdxI)

  Nz = size(F,3)
  N = size(F,1)
  NF = size(F,4)

  zThreads = 30
  threads = (N, N, zThreads)
  zBlocks = cld(N*N*Nz*NF,N*N*zThreads)
  blocks  = (1, 1, zBlocks)
  @. FRhoTrC = 0
  @show blocks,threads
  @cuda blocks=blocks threads=threads DivRhoTrColumnKernel!(FRhoTrC,uC,vC,w,RhoTrC,
    D,dXdxI,N,Nz,NF)
  CUDA.@time @cuda blocks=blocks threads=threads DivRhoTrColumnKernel!(FRhoTrC,uC,vC,
    w,RhoTrC,D,dXdxI,N,Nz,NF)

end

function DivRhoTrColumnKernel!(FRhoTrC,uC,vC,w,RhoTrC,D,dXdxI,N,Nz,NF)

  i = (blockIdx().x-1) * blockDim().x + threadIdx().x
  j = (blockIdx().y-1) * blockDim().y + threadIdx().y
  l = (blockIdx().z-1) * blockDim().z + threadIdx().z
  iz = mod(l,Nz)
  if iz == 0
    iz = Nz
  end
  iF = div(l-iz,Nz) + 1
  
  @inbounds if l <= Nz * NF
    temp1 = -RhoTrC[i,j,iz,iF] * (dXdxI[i,j,1,iz,iF,3,1] * uC[i,j,iz,iF] +
      dXdxI[i,j,1,iz,iF,3,2] * vC[i,j,iz,iF] + dXdxI[i,j,1,iz,iF,3,3] * w[i,j,iz,iF])
    temp2 = -RhoTrC[i,j,iz,iF] * (dXdxI[i,j,2,iz,iF,3,1] * uC[i,j,iz,iF] +
      dXdxI[i,j,2,iz,iF,3,2] * vC[i,j,iz,iF] + dXdxI[i,j,2,iz,iF,3,3] * w[i,j,iz+1,iF])
    CUDA.@atomic FRhoTrC[i,j,iz,iF] += (temp2 - temp1)
    if iz > 1
      CUDA.@atomic FRhoTrC[i,j,iz-1,iF] -= 1/2 * temp1
      CUDA.@atomic FRhoTrC[i,j,iz,iF] += 1/2 * temp2
    end
    if iz < Nz
      CUDA.@atomic FRhoTrC[i,j,iz,iF] -= 1/2 * temp2
      CUDA.@atomic FRhoTrC[i,j,iz+1,iF] += 1/2 * temp1
    end
    tempx = -RhoTrC[i,j,iz,iF] * ((dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * uC[i,j,iz,iF] +
      (dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * vC[i,j,iz,iF] +
      dXdxI[i,j,1,iz,iF,1,3] * w[i,j,iz,iF] + dXdxI[i,j,2,iz,iF,1,3] * w[i,j,iz+1,iF])
    tempy = -RhoTrC[i,j,iz,iF] * ((dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * uC[i,j,iz,iF] +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * vC[i,j,iz,iF] +
      dXdxI[i,j,1,iz,iF,2,3] * w[i,j,iz,iF] + dXdxI[i,j,2,iz,iF,2,3] * w[i,j,iz+1,iF])
    for k = 1 : N
      CUDA.@atomic FRhoTrC[k,j,iz,iF] += D[k,i] * tempx
      CUDA.@atomic FRhoTrC[i,k,iz,iF] += D[k,j] * tempy
    end
  end
  return nothing
end

function DivRhoGradKernelCUDA!(F,cC,RhoC,D,DW,dXdxI,J,N,Nz,NF)

  i = (blockIdx().x-1) * blockDim().x + threadIdx().x
  j = (blockIdx().y-1) * blockDim().y + threadIdx().y
  l = (blockIdx().z-1) * blockDim().z + threadIdx().z
  iz = mod(l,Nz)
  if iz == 0
    iz = Nz
  end
  iF = div(l-iz,Nz) + 1
  
  @inbounds if l <= Nz * NF
    Dxc = 0
    Dyc = 0
    for k = 1 : N
      Dxc = Dxc + D[i,k] * (cC[k,j,iz,iF] / RhoC[k,j,iz,iF])
      Dyc = Dyc + D[j,k] * (cC[i,k,iz,iF] / RhoC[i,k,iz,iF])
    end
    GradDx = ((dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * Dxc +
      (dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * Dyc) / (J[i,j,1,iz,iF] + J[i,j,2,iz,iF])
    GradDy = ((dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * Dxc +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * Dyc) / (J[i,j,1,iz,iF] + J[i,j,2,iz,iF])
    tempx = (dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * GradDx +
      (dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * GradDy
    tempy = (dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * GradDx +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * GradDy
    for k = 1 : N
      CUDA.@atomic F[k,j,iz,iF] = F[k,j,iz,iF] + DW[k,i] * tempx
      CUDA.@atomic F[i,k,iz,iF] = F[i,k,iz,iF] + DW[k,j] * tempy
    end
  end
  return nothing
end

function DivRhoGradCUDA!(F,cC,RhoC,DD,DW,dXdxI,J)
  Nz = size(F,3)
  N = size(F,1)
  NF = size(F,4)
  zThreads = 40
  threads = (N, N, zThreads)
  zBlocks = cld(N*N*Nz*NF,N*N*zThreads) 
  blocks  = (1, 1, zBlocks)
  @. F = 0
  @cuda blocks=blocks threads=threads DivRhoGradKernelCUDA!(F,cC,RhoC,DD,DW,dXdxI,J,N,Nz,NF)
  CUDA.@time @cuda blocks=blocks threads=threads DivRhoGradKernelCUDA!(F,cC,RhoC,DD,DW,dXdxI,J,N,Nz,NF)
end

@kernel function DivRhoGradKernel!(F,@Const(cC),@Const(RhoC),@Const(D),@Const(DW),@Const(dXdxI),@Const(J),N,Nz,NF)


  i, j, iz, iF = @index(Global, NTuple)

  @inbounds begin
    Dxc = 0
    Dyc = 0
    for k = 1 : N
      Dxc = Dxc + D[i,k] * (cC[k,j,iz,iF] / RhoC[k,j,iz,iF])
      Dyc = Dyc + D[j,k] * (cC[i,k,iz,iF] / RhoC[i,k,iz,iF])
    end
    GradDx = ((dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * Dxc +
      (dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * Dyc) / (J[i,j,1,iz,iF] + J[i,j,2,iz,iF])
    GradDy = ((dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * Dxc +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * Dyc) / (J[i,j,1,iz,iF] + J[i,j,2,iz,iF])
    tempx = (dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * GradDx +
      (dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * GradDy
    tempy = (dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * GradDx +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * GradDy
    for k = 1 : N
      @atomic F[k,j,iz,iF] += DW[k,i] * tempx
      @atomic F[i,k,iz,iF] += DW[k,j] * tempy
    end
  end
end

# Creating a wrapper kernel for launching with error checks
function DivRhoGradGPU!(F,cC,RhoC,DD,DW,dXdxI,J)

  Nz = size(F,3)
  N = size(F,1)
  NF = size(F,4)
  zThreads = 40
  groups  = (N, N, zThreads)
  @. F = 0
  backend = KernelAbstractions.get_backend(F)
  kernel! = DivRhoGradKernel!(backend,groups)
  kernel!(F,cC,RhoC,DD,DW,dXdxI,J,N,Nz,NF,ndrange=size(F))
  CUDA.@time kernel!(F,cC,RhoC,DD,DW,dXdxI,J,N,Nz,NF,ndrange=size(F))
end

@kernel function CopyKernel!(F,@Const(c),::Val{BANK}=Val(1)) where BANK

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim)

  # private variable for tile output
  outval = @private eltype(F) 1
  @inbounds outval[1] = -zero(eltype(F))

  for t in 0:ColumnTiles - 1
    Iz = (gz - 1) * ColumnTilesDim + iz
    if Iz <= Nz
      @inbounds cCol[I,J,iz] = c[I,J,Iz,IF]
    else
      @inbounds cCol[I,J,iz] = 0.0
    end
    # wait for all tiles to be loaded
    out = zero(eltype(F))
    out += cCol[I,J,iz]
    outval[1] = out
    @synchronize
    Iz = (gz - 1) * ColumnTilesDim + iz
    if Iz <= Nz
      @inbounds F[I,J,Iz,IF] = outval[1]
    end
  end
end


function Test()
backend = CUDABackend()
#CUDA.allowscalar(false)

KernelAbstractions.synchronize(backend)

Nz = 80
NF = 200
OrdPoly = 4
N = OrdPoly + 1

c = rand!(allocate(backend, Float32, N, N, Nz, NF))
FF = KernelAbstractions.zeros(backend, Float32, N, N, Nz ,NF)

group = (N, N, 40, 1)
kern = CopyKernel!(backend, group)
kern(FF, c, ndrange=size(FF))
KernelAbstractions.synchronize(backend)
CUDA.@time kern(FF, c, ndrange=size(FF))
KernelAbstractions.synchronize(backend)
@show sum(abs.(FF-c))

(DW,DS,) = CGDycore.DerivativeMatrixSingle(OrdPoly)
DW32 = zeros(Float32,OrdPoly+1,OrdPoly+1)
DS32 = zeros(Float32,OrdPoly+1,OrdPoly+1)
@. DW32 = Float32(DW)
@. DS32 = Float32(DS)
F = zeros(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
cC = rand(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
uC = rand(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
vC = rand(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
w = rand(Float32,OrdPoly+1,OrdPoly+1,Nz+1,NF)
RhoC = ones(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
dXdxI = rand(Float32,OrdPoly+1,OrdPoly+1,2,Nz,NF,3,3)
J = ones(Float32,OrdPoly+1,OrdPoly+1,2,Nz,NF)

dDW = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1)
dDS = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1)
copyto!(dDW,DW32)
copyto!(dDS,DS32)
dF = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,Nz,NF) 
dcC = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
duC = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
dvC = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
dw = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,Nz+1,NF)
dRhoC = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
ddXdxI = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,2,Nz,NF,3,3)
dJ = KernelAbstractions.zeros(backend,Float32,OrdPoly+1,OrdPoly+1,2,Nz,NF)
copyto!(dF,F)
copyto!(dcC,cC)
copyto!(duC,uC)
copyto!(dvC,vC)
copyto!(dw,w)
copyto!(dRhoC,RhoC)
copyto!(ddXdxI,dXdxI)
copyto!(dJ,J)

DivRhoGradGPU!(dF,dcC,dRhoC,dDS,dDW,ddXdxI,dJ)
KernelAbstractions.synchronize(backend)
copyto!(F,dF)
@show "Kern",sum(abs.(F));

N = OrdPoly + 1
@. F = 0
for iF = 1 : NF
  @views CGDycore.DivRhoGradE!(F[:,:,:,iF],cC[:,:,:,iF],RhoC[:,:,:,iF],DS32,DW32,
    dXdxI[:,:,:,:,iF,:,:],J[:,:,:,:,iF])
end
@time for iF = 1 : NF
  @views CGDycore.DivRhoGradE!(F[:,:,:,iF],cC[:,:,:,iF],RhoC[:,:,:,iF],DS32,DW32,
    dXdxI[:,:,:,:,iF,:,:],J[:,:,:,:,iF])
end
@show "CPU",sum(abs.(F))
DivRhoGradCUDA!(dF,dcC,dRhoC,dDS,dDW,ddXdxI,dJ)
copyto!(F,dF)
@show "CUDA",sum(abs.(F));
stop

@. F = 0
for iF = 1 : NF
  @views DivRhoTrColumnE1!(F[:,:,:,iF],uC[:,:,:,iF],vC[:,:,:,iF],w[:,:,:,iF],
    cC[:,:,:,iF],DS32,dXdxI[:,:,:,:,iF,:,:])
end
@time for iF = 1 : NF
  @views DivRhoTrColumnE1!(F[:,:,:,iF],uC[:,:,:,iF],vC[:,:,:,iF],w[:,:,:,iF],
    cC[:,:,:,iF],DS32,dXdxI[:,:,:,:,iF,:,:])
end
@show "CPU",sum(abs.(F))
DivRhoTrColumnE2!(dF,duC,dvC,dw,dcC,dDS,ddXdxI)
copyto!(F,dF)
@show "GPU",sum(abs.(F));
end
