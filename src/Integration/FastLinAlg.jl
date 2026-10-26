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

# y += sum_i a_i x_i with Kahan compensation c (one array of the size of y).
# In Float32 the increments of large components (rho theta ~ 300) can be
# smaller than one ulp and are otherwise lost.
function AXPYKahan!(y,c,x,a,Global)
  n = length(a)
  if n == 0
    return
  end
  backend = get_backend(y)
  ly = length(y)
  aS = SVector{length(a)}(a)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  yV = reshape(y,length(y))
  cV = reshape(c,length(y))
  xV = reshape(x,length(y),n)
  group = (NumberThreadGPU)
  ndrange = (ly)
  KAXPYKahanKernel! = AXPYKahanKernel!(backend,group)
  KAXPYKahanKernel!(yV,cV,xV,aS,Val(n),ndrange=ndrange)
end

@kernel inbounds = true function AXPYKahanKernel!(y, c, @Const(x), @Const(a), ::Val{n}) where n
  I, = @index(Global, NTuple)

  N = @uniform @ndrange()[1]

  if I <= N
    sum_val = eltype(y)(0)
    @unroll for i in 1:n
      sum_val += a[i] * x[I,i]
    end
    yv = sum_val - c[I]
    t = y[I] + yv
    c[I] = (t - y[I]) - yv
    y[I] = t
  end
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

