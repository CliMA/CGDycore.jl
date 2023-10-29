function cart2sphere(x,y,z)
FT = eltype(x)
@show FT
r=sqrt(x^2+y^2+z^2);
phi=asin(z/r);

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=FT(0)
if abs(abs(phi)-FT(pi/2))>FT(1.e-14)
  lam=atan(y,x); # TODO: check translation with Oswald
  if lam<FT(0.0)
    lam=lam+FT(2*pi);
  end
end
return (lam,phi,r)
end

function cart2sphereDeg(x,y,z)
FT = eltype(x)
r=sqrt(x^2+y^2+z^2);
phi=asind(z/r);

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=-FT(180);
if abs(abs(phi)-FT(90.0))>FT(1.e-14)
  lam=atand(y,x); # TODO: check translation with Oswald
end
return (float(lam),float(phi),r)
end


