function T=TFromThetaRho(Th,Rho,Param)
T=(Param.Rd*Th/Param.p0^Param.kappa).^(1/(1-Param.kappa))./(Param.Rd*Rho);
end