using MPI
using Base
using CUDA
using AMDGPU
using KernelAbstractions
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras



@kernel function FluxUpdateKernel!(Flux,@Const(u),@Const(v),@Const(w),@Const(c))

  i,j,k = @index(Global, NTuple)

  M = @uniform @ndrange()[1]
  N = @uniform @ndrange()[2]
  L = @uniform @ndrange()[3]

  if 1 < i && i < M && j > 1 &&  j < N && k > 1 && k < L
    @inbounds FluxW = eltype(Flux)(0.5) * (u[i-1,j,k] - abs(u[i-1,j,k])) * c[i,j,k] 
    @inbounds FluxE = eltype(Flux)(0.5) * (u[i,j,k] + abs(u[i,j,k])) * c[i,j,k] 
    @inbounds FluxS = eltype(Flux)(0.5) * (v[i,j-1,k] - abs(v[i,j-1,k])) * c[i,j,k] 
    @inbounds FluxN = eltype(Flux)(0.5) * (v[i,j,k] + abs(v[i,j,k])) * c[i,j,k] 
    @inbounds FluxB = eltype(Flux)(0.5) * (w[i,j,k-1] - abs(w[i,j,k-1])) * c[i,j,k] 
    @inbounds FluxT = eltype(Flux)(0.5) * (w[i,j,k] + abs(w[i,j,k])) * c[i,j,k] 
#   @inbounds @atomic Flux[i-1,j,k] += -FluxW
#   @inbounds @atomic Flux[i+1,j,k] += FluxE
#   @inbounds @atomic Flux[i,j-1,k] += -FluxS
#   @inbounds @atomic Flux[i,j+1,k] += FluxN
#   @inbounds @atomic Flux[i,j,k-1] += -FluxB
#   @inbounds @atomic Flux[i,j,k+1] += FluxT
    ind = i + j + k
    @inbounds @atomic Flux[ind] += (FluxW - FluxE + FluxS - FluxN + FluxB - FluxT)
  end  
end

backend = CUDABackend()
FT = Float32
M = 100
N = 101
L = 102

group = (9, 9 ,9)
ndrange = (M, N , L)

u = KernelAbstractions.ones(backend,FT,M+1,N,L)
v = KernelAbstractions.ones(backend,FT,M,N+1,L)
w = KernelAbstractions.ones(backend,FT,M,N,L+1)
c = KernelAbstractions.ones(backend,FT,M,N,L)
Flux = KernelAbstractions.zeros(backend,FT,M+N,L)

KFluxUpdateKernel! = FluxUpdateKernel!(backend,group)

KFluxUpdateKernel!(Flux,u,v,w,c,ndrange=ndrange)
@time for iter = 1 : 100
  KFluxUpdateKernel!(Flux,u,v,w,c,ndrange=ndrange)
end  
