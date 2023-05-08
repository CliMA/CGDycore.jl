using CGDycore
using Metal
using LinearAlgebra

a = ones(Float32,30,30)
b = ones(Float32,30,30)

ad = MtlArray(a)
bd = MtlArray(b)
cd = similar(ad)

@metal threads=length(cd) CGDycore.vadd(ad, bd, cd)

dims = tuple(30,30)
# Create a Metal array with a default storage mode of shared (both CPU and GPU get access)
arr_mtl = MtlArray{Float32}(undef, dims)
brr_mtl = MtlArray{Float32}(undef, dims)
crr_mtl = MtlArray{Float32}(undef, dims)
# Unsafe wrap the contents of the Metal array with a CPU array
arr_cpu = unsafe_wrap(Array{Float32}, arr_mtl, dims)
brr_cpu = unsafe_wrap(Array{Float32}, brr_mtl, dims)
crr_cpu = unsafe_wrap(Array{Float32}, crr_mtl, dims)
@. arr_cpu = 4
@. brr_cpu = 4
brr_cpu[4,4] = 8
@metal threads=length(crr_mtl) CGDycore.vmult(arr_mtl, brr_mtl, crr_mtl)
@show crr_cpu[4,4]
@time mul!(cd,ad,bd)
@time mul!(cd,ad,bd)
@time mul!(cd,ad,bd)
c=Array(cd)
@time mul!(c,a,b)
@time mul!(c,a,b)
@time mul!(c,a,b)
@show c[1,1]
