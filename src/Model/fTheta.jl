function Th = fTheta(x,Param)
switch lower(Param.ProfTheta)
  case 'solidbody'
    [Lon,Lat,R]=cart2sphere(x(1),x(2),x(3));
    z=max(R-Param.RadEarth,0);
    NBr=Param.NBr;
    Grav=Param.Grav;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    Th=Th0*exp(z*S);
  case 'hyperdiffcart'
    if abs(x(1)-3000.*1.e3)<=1000*1.e3 && abs(x(2))<=1000*1.e3
      Th=100.;
    else
      Th=0;
    end
  case 'hyperdiff'
    [Lon,Lat,R]=cart2sphere(x(1),x(2),x(3));
    GreatCircleR=acos(sin(Param.PertLat)*sin(Lat)...
      +cos(Param.PertLat)*cos(Lat)*cos(Lon-Param.PertLon));
    GreatCircleR=GreatCircleR/Param.PertExpR;
    % Apply perturbation in zonal velocity
    if GreatCircleR<6.0
      Th=sqrt(2*Param.ValDiff);
    else
      Th=0.0;
    end
  case 'barowavecart'
    eta=EtaFromZ(x(1),x(2),x(3),Param);
    p=Param.p0*eta;
    T=TBaroWave(x(1),x(2),eta,Param);
    Th=T*(Param.p0/p)^(Param.Rd/Param.Cpd);
  case 'barowavesphere'
    [Lon,Lat,R]=cart2sphere(x(1),x(2),x(3));
    Z=max(R-Param.RadEarth,0);
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Param.Rd*T0/Param.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z)...
      +ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0)...
      +ConstB*Z*exp(-ScaledZ*ScaledZ);
    IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
    if Param.Deep
      RRatio= R/Param.EarthRadius;
    else
      RRatio = 1.0;
    end
    InteriorTerm=(RRatio*cos(Lat))^Param.K...
      -Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0);
    Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);
    Pressure=Param.p0*exp(-Param.Grav/Param.Rd...
      *(IntTau1-IntTau2*InteriorTerm));
    Th=Temperature*(Param.p0/Pressure)^(Param.Rd/Param.Cpd);
  case 'baldaufsphere'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    r=r-Param.RadEarth;
    p=Param.p0*exp(-Param.Grav*r/(Param.Rd*Param.T0));
    Rho=p/(Param.Rd*Param.T0);
    T=p/(Rho*Param.Rd);
    d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    T=T+Param.DeltaT*exp(-Param.ExpDist*d)*sin(pi*r/Param.H);
    Th=T*(Param.p0/p)^(Param.Rd/Param.Cpd);
  case 'baldaufcart'
    delta=Param.Grav/(Param.Rd*Param.T0);
    p=Param.p0*exp(-delta*x(3));
    
    dT=Param.DeltaT*exp(-(x(1)-Param.xc)^2/Param.d^2)*sin(pi*x(3)/Param.H);
    
    TLoc=Param.T0+exp(delta/2*x(3))*dT;
    Th=TLoc*(Param.p0/p)^(Param.Rd/Param.Cpd);  
  case 'warmbubble2d'
    Th0=Param.Th0;
    DeltaTh=Param.DeltaTh;
    xC0=Param.xC0;
    zC0=Param.zC0;
    rC0=Param.rC0;
    x3=x(3);
    x1=x(1);
    rr=sqrt((x1-xC0)^2+(x3-zC0)^2);
    ThLoc=Th0;
    if rr<rC0
      ThLoc=ThLoc+DeltaTh*cos(0.5*pi*rr/rC0)^2;
    end
    Th=ThLoc;
  case 'densitycurrent'
    T0=Param.T0;
    DeltaT=Param.DeltaT;
    xC0=Param.xC0;
    zC0=Param.zC0;
    xrC0=Param.xrC0;
    zrC0=Param.zrC0;
    x3=x(3);
    x1=x(1);
    pLoc=Param.p0*(1-Param.kappa*Param.Grav*x(3)...
      /(Param.Rd*T0))^(Param.Cpd/Param.Rd);
    Rad=sqrt(((x1-xC0)/xrC0)^2+((x3-zC0)/zrC0)^2);
    Th=T0;
    if Rad<1.0d0
      Th=Th+DeltaT*(cos(pi*Rad)+1.0)/2.0*(pLoc/Param.p0)^(-Param.kappa);
    end
  case 'gravityhill'
    z=x(3);
    NBr=Param.NBr;
    Grav=Param.Grav;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    Th=Th0*exp(z*S);
  case 'inertiagravitywave'
    z=x(3);
    NBr=Param.NBr;
    Grav=Param.Grav;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    ThB=Th0*exp(z*S);
    Th=ThB+Param.DeltaTh*sin(pi*x(3)/Param.H)/(1+(x(1)-Param.xC)^2/Param.a^2);
  case 'galewsky'
    Th=1;
  case 'rossbyhaurwitz'
    Grav=Param.Grav;
    Omega=Param.Omega;
    H06=8000.0;
    omega6=7.8480e-6;
    K6=7.8480e-6;
    R6=4.0;
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    A=0.5*omega6*(2.0*Omega+omega6)*cos(lat)*cos(lat)...
      +0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat)...
      +(2.0*R6*R6-R6-2.0)-2.0*R6*R6*(cos(lat))^(-2));
    B=2.0*(Omega+omega6)*K6/(R6+1.0)/(R6+2.0)...
      *(cos(lat))^R6*((R6*R6+2.0*R6+2.0)...
      -(R6+1.0)*(R6+1.0)*cos(lat)*cos(lat));
    C=0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat)-(R6+2.0));
    Th=(Grav*H06+r*r*(A+B*cos(R6*lon)+C*cos(2.0*R6*lon)))/Grav;
  case 'hill'
    lon0=3.0e0/2.0e0*pi;
    lat0=pi/6.0e0;
    H05=5960.0e0;
    hS=2000.0e0;
    RadiusC=pi/9.0e0;
    UMax=20.0e0;
    rotation_angle=0;
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    [rot_lon,rot_lat]=Rotate(lon,lat,rotation_angle);
    [rot_lon0,rot_lat0]=Rotate(lon0,lat0,rotation_angle);
    r=sqrt(min(RadiusC*RadiusC, ...
      (rot_lon-rot_lon0)*(rot_lon-rot_lon0)+(rot_lat-rot_lat0)*(rot_lat-rot_lat0)));
    HeightLoc=hS*(1-r/RadiusC);
    Th=(Param.Grav*H05-(Param.RadEarth*Param.Omega*UMax+0.5*UMax*UMax)*sin(rot_lat)*sin(rot_lat))...
      /Param.Grav-HeightLoc;
  case 'spherical'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    if abs(d)<=0.8
      Th=100*cos(pi*d/0.8/2)^2+100;
    else
      Th=100.0;
    end
  case 'cosinebell'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    r1=Param.RadEarth*acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    R3=Param.RadEarth/3.0d0;
    if r1<=R3
      Th=1000.0e0/2.0e0*(1.0d0+cos(pi*r1/R3));
    else
      Th=0.0d0;
    end
  case 'sphericalsmooth'
    Th=exp(-((-x*normal-1.0)/0.1)^2);
    %Th=exp(-((-x(2)-1.0)/0.4)^2)
  case 'quad'
    if abs(x(1)-xM(1))<xH(1) && abs(x(2)-xM(2))<xH(2)
      Th=1;
    else
      Th=0;
    end
  case 'cosine'
    d=sqrt((x(1)-xM(1))^2+(x(2)-xM(2))^2);
    if d<=rH
      Th=cos(pi*d/rH/2)^2;
    else
      Th=0.0;
    end
  case 'constant'
    Th=1;
  case 'linear'
    Th=x(1)+1;
  case 'bickley'
    Th=Param.RhoTheta;
end
end

function intG=integrandG(tau,RadEarth)
global Omega uM lat0G lat1G eN
f=2.0*Omega*sin(tau);
if (tau<=lat0G) || (tau>=lat1G)
  uStart=0.0;
else
  uStart=uM/eN*exp(1.0/((tau-lat0G)*(tau-lat1G)));
end
if abs(tau)<0.5*pi
  intG=(RadEarth*f+uStart*tan(tau))*uStart;
else
  intG=0.0;
end
end

function [rot_lon,rot_lat]=Rotate(lon,lat,rotation_angle)
if abs(rotation_angle)<1.0e-8
  rot_lon = lon;
  rot_lat = lat;
else
  [rot_lon,rot_lat]=regrot(lon,lat,0.0e0,-0.5e0*pi+rotation_angle);
end
end

function [pxrot,pyrot]=regrot(pxreg,pyreg,pxcen,pycen)
  
%----------------------------------------------------------------------
%
%*    conversion between regular and rotated spherical coordinates.
%*
%*    pxreg     longitudes of the regular coordinates
%*    pyreg     latitudes of the regular coordinates
%*    pxrot     longitudes of the rotated coordinates
%*    pyrot     latitudes of the rotated coordinates
%*              all coordinates given in degrees n (negative for s)
%*              and degrees e (negative values for w)
%*    pxcen     regular longitude of the south pole of the rotated grid
%*    pycen     regular latitude of the south pole of the rotated grid
%*
%*    kcall=-1: find regular as functions of rotated coordinates.
%*    kcall= 1: find rotated as functions of regular coordinates.
%
%-----------------------------------------------------------------------
%

zpih = pi*0.5d0;
zsycen = sin((pycen+zpih));
zcycen = cos((pycen+zpih));

zxmxc  = pxreg - pxcen;
zsxmxc = sin(zxmxc);
zcxmxc = cos(zxmxc);
zsyreg = sin(pyreg);
zcyreg = cos(pyreg);
zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc;
zsyrot = max(zsyrot,-1.e0);
zsyrot = min(zsyrot,+1.e0);

pyrot = asin(zsyrot);

zcyrot = cos(pyrot);
zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot;
zcxrot = max(zcxrot,-1.e0);
zcxrot = min(zcxrot,+1.e0);
zsxrot = zcyreg*zsxmxc/zcyrot;

pxrot = acos(zcxrot);

if zsxrot<0.0
  pxrot = -pxrot;
end
end

function [pures,pvres]=turnwi(puarg,pvarg, ...
                  pxreg,pyreg,pxrot,pyrot, ...
                  pxcen,pycen)
%
%-----------------------------------------------------------------------
%
%*    turn horizontal velocity components between regular and
%*    rotated spherical coordinates.
%
%*    puarg : input u components
%*    pvarg : input v components
%*    pures : output u components
%*    pvres : output v components
%*    pa    : transformation coefficients
%*    pb    :    -"-
%*    pc    :    -"-
%*    pd    :    -"-
%*    pxreg : regular longitudes
%*    pyreg : regular latitudes
%*    pxrot : rotated longitudes
%*    pyrot : rotated latitudes
%*    kxdim              : dimension in the x (longitude) direction
%*    kydim              : dimension in the y (latitude) direction
%*    kx                 : number of gridpoints in the x direction
%*    ky                 : number of gridpoints in the y direction
%*    pxcen              : regular longitude of the south pole of the
%*                         transformed grid
%*    pycen              : regular latitude of the south pole of the
%*                         transformed grid
%*
%*    kcall < 0          : find wind components in regular coordinates
%*                         from wind components in rotated coordinates
%*    kcall > 0          : find wind components in rotated coordinates
%*                         from wind components in regular coordinates
%*    note that all coordinates are given in degrees n and degrees e.
%*       (negative values for s and w)
%
!-----------------------------------------------------------------------

     zpih = pi*0.5e0;
     zsyc = sin(pycen+zpih);
     zcyc = cos(pycen+zpih);
     
     zsxreg = sin(pxreg);
     zcxreg = cos(pxreg);
     zsyreg = sin(pyreg);
     zcyreg = cos(pyreg);
     
     zxmxc  = pxreg - pxcen;
     zsxmxc = sin(zxmxc);
     zcxmxc = cos(zxmxc);
     
     zsxrot = sin(pxrot);
     zcxrot = cos(pxrot);
     zsyrot = sin(pyrot);
     zcyrot = cos(pyrot);
    
     pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot;
     pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot - ...
          zcxmxc*zsxrot*zsyrot;
     pc =-zsyc*zsxrot/zcyreg;
     pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg;
     !
     pures = pa*puarg + pb*pvarg;
     pvres = pc*puarg + pd*pvarg;
end




