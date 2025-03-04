function TestShallowWaterDG2(Example)
close all
% Number of grid cells per panel 
n=10;
% DG order and DG integration order
OrdPolyX=4;
%Number of variables
nV=3;
hPos=1;
uPos=2;
vPos=3;
wPos=4;

%Physical parameters
Param.Grav=9.81d0;
Param.Cpd=1004.0d0;
Param.Cvd=717.0d0;
Param.Rd=Param.Cpd-Param.Cvd;
Param.p0=1.0d5;
Param.Cpv=1885.0d0;
Param.Gamma=Param.Cpd/Param.Cvd;
Param.kappa=Param.Rd/Param.Cpd;
cS=340;
cS=sqrt(Param.Grav*1.e5);
lambda=cS;
Param.cS=cS;
Param.lambda=Param.cS;

U=.0;
beta=1;
NBr=1.e-2;
Param.b0=0.01;
Param.A=5000;
Param.H=10000;
Param.xC=0;
Param.Th0=300;
Param.DeltaTh=2;


%Output
Subs=4;
fig=1;
OutputStart=true;
OutputEnd=true;
% Gitter
nS=32;
nRad=4;

Param.Grav=9.81;
Param.Omega=2*pi/(24*3600);
Param.RadEarth=6.37122d+6;
Grid=CubedGrid(n,@OrientFaceSphere,Param);

Param.OrdPolyX=OrdPolyX;
Param.nV=nV;
Param.hPos=hPos;
Param.uPos=uPos;
Param.vPos=vPos;
Param.wPos=wPos;

Param.Profh=Example;
Param.ProfuVel=Example;

Param.IntMethod='RungeKutta';
RKMethod='RK4';
Param.RK=RungeKuttaMethod(RKMethod);
Param.GridType='Radial';
Param.OutputEnd=true;
[V,Param]=InitDG2(@JacobiSphere2Cart,@FcnNonLin,Grid,Param);

end



      


