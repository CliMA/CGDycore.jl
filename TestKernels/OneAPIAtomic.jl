using KernelAbstractions
using Atomix
using Adapt
using oneAPI

@kernel function mykernel(x, y)
    i = @index(Global)
    Atomix.@atomic x[1] += y[i]
end

backend = oneAPIBackend()

x_ = zeros(Float32, 1)
y_ = rand(Float32, 100)
x = Adapt.adapt(backend, x_)
y = Adapt.adapt(backend, y_)

mykernel(backend)(x, y, ndrange=length(y))

@info "" sum(x) sum(y)
