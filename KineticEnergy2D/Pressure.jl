function Pressure!(p,RhoTh,Phys)
  p0 = Phys.p0
  Rd = Phys.Rd
  kappa = Phys.kappa
  @inbounds for i in eachindex(p)  
    p[i] = p0 * fast_pow(Rd * RhoTh[i] / p0, 1.0 / (1.0 - kappa));
  end 
end 

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
fast_pow(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))
