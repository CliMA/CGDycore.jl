function fpBGrd(x,Param)
if lowercase(Param.ProfpBGrd)== "baldaufcart"
    delta=Param.Grav/(Param.Rd*Param.T0);

    p=Param.p0*exp(-delta*x[3]);
else
    p=0;
end
return p
end
