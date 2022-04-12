function cart2sphere(x,y,z)
r=sqrt(x^2+y^2+z^2);
phi=asin(z/r);

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=0;
if abs(abs(phi)-pi/2)>1.e-14
  lam=atan(y,x); # TODO: check translation with Oswald
  if lam<0.0
    lam=lam+2*pi;
  end
end
return (float(lam),float(phi),r)
end

function cart2sphereDeg(x,y,z)
r=sqrt(x^2+y^2+z^2);
phi=asind(z/r);

# ϕ = atan(z, hypot(y, x))
# if abs(ϕ) == 90
#     λ = zero(ϕ)
# else
#     λ = atan(y, x)
# end

lam=-180;
if abs(abs(phi)-90.0)>1.e-14
  lam=atand(y,x); # TODO: check translation with Oswald
end
return (float(lam),float(phi),r)
end


