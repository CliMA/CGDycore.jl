import KernelAbstractions
using Metal
using MetalKernels


backend = MtlBackend()
function vadd(a, b, c)
  i = thread_position_in_grid_1d()
  c[i] = a[i] + b[i]
  return
end

@kernel function mul2_kernel(A)
  I = @index(Global)
  A[I] = 2 * A[I]
end

A = MtlArray(ones(Float32,1024, 1024))
mul2_kernel(backend, 64)(A, ndrange=size(A))
synchronize(backend)
all(A .== 2.0)

dev = CPU()
B = ones(1024, 1024)
ev = mul2_kernel(dev, 64)(B, ndrange=size(B))
synchronize()
@show B[1,1]
all(B .== 2.0)

