function Param=PhysParameters()
% Physical parameters
Param.Grav=9.81; 
Param.cS=360;
Param.Grav=9.81d0;
Param.Cpd=1004.0d0;
Param.Cvd=717.0d0;
Param.Rd=Param.Cpd-Param.Cvd;
Param.p0=1.0d5;
Param.Cpv=1885.0d0;
Param.Gamma=Param.Cpd/Param.Cvd;
Param.kappa=Param.Rd/Param.Cpd;

% Sphere
Param.Omega=2*pi/(24*3600);
Param.RadEarth=6.37122d+6;
end