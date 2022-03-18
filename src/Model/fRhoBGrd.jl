function fRhoBGrd(x,Param)
if lower(Param.ProfRhoBGrd)== "baldaufcart"
    delta=Param.Grav/(Param.Rd*Param.T0);
    p=Param.p0*exp(-delta*x[3]);
    TLoc=Param.T0;
    Rho=p/(Param.Rd*TLoc);
else
    Rho=0;
end
return Rho
end
