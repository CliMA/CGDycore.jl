function hS(x,Param)
f=0;
if Param.hS == "AgnesiAnnulus"
    f=Param.h/(((x-Param.phi0)/Param.aphi)^2+1);
elseif Param.hS == "AgnesiCart"
    f=Param.hC/(((x-Param.x0C)/Param.aC)^2+1);
else
    f=0;
end
return f
end