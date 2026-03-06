function AXPY!(y,x,a,Global)
  n = length(a)
  if n == 0 
    return
  end    
  backend = get_backend(y)
  ly = length(y)
  aS = SVector{length(a)}(a)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  yV = reshape(y,length(y))
  xV = reshape(x,length(y),n)
  group = (NumberThreadGPU)
  ndrange = (ly)
  KAXPYKernel! = AXPYKernel!(backend,group)
  KAXPYKernel!(yV,xV,aS,Val(n),ndrange=ndrange)
end

function AXPY!(y,x,dt,a,Global)
  n = length(a)
  if n == 0 
    return
  end    
  backend = get_backend(y)
  ly = length(y)
  aS = SVector{length(a)}(dt*a)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  yV = reshape(y,length(y))
  xV = reshape(x,length(y),n)
  group = (NumberThreadGPU)
  ndrange = (ly)
  KAXPYKernel! = AXPYKernel!(backend,group)
  KAXPYKernel!(yV,xV,aS,Val(n),ndrange=ndrange)
end

@kernel inbounds = true function AXPYKernel!(y, @Const(x), @Const(a), ::Val{n}) where n
  I, = @index(Global, NTuple)

  N = @uniform @ndrange()[1]
    
  if I <= N
    sum_val = eltype(y)(0)
    @unroll for i in 1:n
      sum_val += a[i] * x[I,i]
    end
    y[I] += sum_val
  end
end

