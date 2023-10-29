function sphere2cart(lam,phi,r)
x=cos(lam)*cos(phi)*r;
y=sin(lam)*cos(phi)*r;
z=sin(phi)*r;
return [x;y;z];
end
function sphereDeg2cart(lam,phi,r)
x=cosd(lam)*cosd(phi)*r;
y=sind(lam)*cosd(phi)*r;
z=sind(phi)*r;
return [x;y;z];
end


