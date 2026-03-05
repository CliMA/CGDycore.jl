function AXPY!(y,x,a,Global)
  n = length(a)
  if n == 0 
    return
  end    
  backend = get_backend(y)
  M = size(y,1)
  Nz = size(y,2)
  ND = size(y,3)
  NV = size(y,4)
  aS = SVector{length(a)}(a)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,ND,NV)
  KAXPYKernel! = AXPYKernel!(backend,group)
  KAXPYKernel!(y,x,aS,Val(n),ndrange=ndrange)
end

function AXPY!(y,x,dt,a,Global)
  n = length(a)
  if n == 0 
    return
  end    
  backend = get_backend(y)
  M = size(y,1)
  Nz = size(y,2)
  ND = size(y,3)
  NV = size(y,4)
  aS = SVector{length(a)}(dt*a)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG,1)
  ndrange = (M,Nz,ND,NV)
  KAXPYKernel! = AXPYKernel!(backend,group)
  KAXPYKernel!(y,x,aS,Val(n),ndrange=ndrange)
end

@kernel inbounds = true function AXPYKernel!(y, @Const(x), @Const(a), ::Val{n}) where n
  K, IZ, ID, IV = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]
  NV = @uniform @ndrange()[4]
    
  if ID <= ND && IV <= NV
    sum_val = a[1] * x[K,IZ,ID,IV,1]
    @unroll for i in 2:n
      sum_val += a[i] * x[K,IZ,ID,IV,i]
    end
    y[K,IZ,ID,IV] += sum_val
  end
end

