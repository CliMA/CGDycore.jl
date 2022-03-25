function sphere2cart(lam,phi,r)
x=cos(lam)*cos(phi)*r;
y=sin(lam)*cos(phi)*r;
z=sin(phi)*r;
return [x;y;z];
end


