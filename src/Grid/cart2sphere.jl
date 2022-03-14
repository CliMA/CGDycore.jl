function [lam,phi,r]=cart2sphere(x,y,z)
r=sqrt(x^2+y^2+z^2);
phi=asin(z/r);
lam=0;
if abs(abs(phi)-pi/2)>1.e-14
  lam=atan2(y,x);
  if lam<0.0
    lam=lam+2*pi;
  end
end
end


