function p=pFromTheta(Th,Param)
p=(Param.Rd*Th/Param.p0^Param.kappa).^(1/(1-Param.kappa));
end