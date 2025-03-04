function p=PresSh(V,Param)
Grav=Param.Grav;
hPos=Param.hPos;
p=0.5*Grav*V(hPos,:,:).^2;
end