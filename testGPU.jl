using CGDycore
using Metal

a = ones(Float32,30,30)
b = ones(Float32,30,30)

ad = MtlArray(a)
bd = MtlArray(b)
cd = similar(ad)

@metal threads=length(cd) vadd(ad, bd, cd)

