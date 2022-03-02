function [uS,vS] = fVel(x,Param)
switch lower(Param.ProfVel)
  case 'solidbody'
    uS=0;
    vS=0;
  case 'hyperdiffcart'
    if abs(x(1)-3000.*1.e3)<=1000*1.e3 && abs(x(2))<=1000*1.e3
      uS=sqrt(5000);
      vS=sqrt(5000);
    else
      uS=0;
      vS=0;
    end
  case 'hyperdiff'
    [Lon,Lat,R]=cart2sphere(x(1),x(2),x(3));
    GreatCircleR=acos(sin(Param.PertLat)*sin(Lat)...
      +cos(Param.PertLat)*cos(Lat)*cos(Lon-Param.PertLon));
    GreatCircleR=GreatCircleR/Param.PertExpR;
    % Apply perturbation in zonal velocity
    if GreatCircleR<6.0
      uS=sqrt(Param.ValDiff);
      vS=sqrt(Param.ValDiff);
    else
      uS=0.0;
      vS=0.0;
    end
  case 'barowavesphere'
    [Lon,Lat,R]=cart2sphere(x(1),x(2),x(3));
    Z=max(R-Param.RadEarth,0);
%     if Z>27500 && abs(Lat-pi/4)<=.1
%       aa=3;
%     end
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Param.Rd*T0/Param.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z)...
    +ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
    if Param.Deep
      RRatio= R/Param.EarthRadius;
    else
      RRatio = 1.0;
    end
    InteriorTerm=(RRatio*cos(Lat))^Param.K...
      -Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0);
    Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);
    
    InteriorTermU=(RRatio*cos(Lat))^(Param.K-1.0)...
      - (RRatio*cos(Lat))^(Param.K+1.0);
    BigU=Param.Grav/Param.RadEarth*Param.K...
      *IntTau2*InteriorTermU*Temperature;
    if Param.Deep
      RCosLat=R*cos(Lat);
    else
      RCosLat=Param.RadEarth*cos(Lat);
    end
    OmegaRCosLat =Param.Omega*RCosLat;
    
    % 		if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
    % 			_EXCEPTIONT("Negative discriminant detected.");
    % 		}
    
    uS=-OmegaRCosLat+sqrt(OmegaRCosLat*OmegaRCosLat+RCosLat*BigU);
    vS=0;
    % Exponential perturbation
    GreatCircleR=acos(sin(Param.PertLat)*sin(Lat)...
      +cos(Param.PertLat)*cos(Lat)*cos(Lon-Param.PertLon));
    
    GreatCircleR=GreatCircleR/Param.PertExpR;
    
    % Tapered perturbation with height
    if Z < Param.PertZ
      PertTaper=1.0-3.0*Z*Z/(Param.PertZ * Param.PertZ)...
        + 2.0*Z*Z*Z/(Param.PertZ*Param.PertZ*Param.PertZ);
    else
      PertTaper=0.0;
    end
    
    % Apply perturbation in zonal velocity
    if GreatCircleR<1.0
      uSPert=Param.Up*PertTaper*exp(-GreatCircleR*GreatCircleR);
    else
      uSPert=0.0;
    end
    uS=uS+uSPert;
    
  case 'barowavecart'
    eta=EtaFromZ(x(1),x(2),x(3),Param);
    uS=-Param.u0*sin(pi*x(2)/Param.Ly)^2*log(eta)...
      *exp(-(log(eta)/Param.b)^2);
    if x(2)<0
      uS=-uS;
    else
      uS=uS+Param.uP*exp(-((x(1)-Param.xC)^2 ...
        +(x(2)-Param.yC)^2)/Param.Lp^2);
    end
    vS=0;
  case 'baldauf'
    uS=0;
    vS=0;
  case 'barowavesphere'
%     # parameters
 
  case 'gravityhill'
    uS=Param.uMax;
    vS=0;
  case 'sphericalharmonics'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    [uS,vS]=harmonicVec(Param.lHar,Param.mHar,lat,lon,r);
  case 'galewsky'
    uM=80.0;
    lat0G=pi/7.0;
    lat1G=pi/2.0-lat0G;
    eN=exp(-4.0/(lat1G-lat0G)^2.0);
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    if (lat<=lat0G) || (lat>=lat1G)
      uS=0;
    else
      uS=uM/eN*exp(1.0/((lat-lat0G)*(lat-lat1G)));
    end
    vS=0;
  case 'rossbyhaurwitz'
    omega6=7.8480e-6; 
    K6=7.8480e-6;
    R6=4.0;
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    uS=r*omega6*cos(lat)...
       +r*K6*(cos(lat))^(R6-1.0)...
       *(R6*sin(lat)*sin(lat)-cos(lat)*cos(lat))*cos(R6*lon);
     vS=-r*K6*R6*(cos(lat))^(R6-1.0)*sin(lat)*sin(R6*lon);
  case 'hill'
    uMax=20;
    rotation_angle=0;
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    [rot_lon,rot_lat]=Rotate(lon,lat,rotation_angle);
    u_lat=uMax*cos(rot_lat);
    v_lat=0.0e0;
    if abs(rotation_angle)<1.0E-8
      uS = u_lat;
      vS = v_lat;
    else
      % rotate wind components
      [uS,vS]= turnwi(u_lat,v_lat,lon,lat,rot_lon,rot_lat,0.0d0,-0.5e0*pi+rotation_angle);
      if abs(uS)<1.0E-10
        uS=0.0e0;
      end
    end
  case 'bickley'
    U1 = cosh(x(2))^(-2);

    % ?? = exp(-(x2 + Param.l / 10)^2 / 2p.l^2) * cos(p.k * x1) * cos(p.k * x2)
    % Vortical velocity fields (u??, u??) = (-?²??, ?¹??)
    gaussian = exp(-(x(2) + Param.l / 10)^2 / (2*Param.l^2));
    u1 = gaussian * (x(2) + Param.l / 10) / Param.l^2 * cos(Param.k * x(1)) * cos(Param.k * x(2));
    u1 = u1+Param.k * gaussian * cos(Param.k * x(1)) * sin(Param.k * x(2));
    u2 = -Param.k * gaussian * sin(Param.k * x(1)) * cos(Param.k * x(2));
    uS=U1+Param.Eps*u1;
    vS=Param.Eps*u2;
  case 'cosinebell' 
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    alpha0=pi/2;
    u0=2*pi*Param.RadEarth/(86400*12);
    uS=u0*(cos(alpha0)*cos(lat)+sin(alpha0)*cos(lon)*sin(lat));
    vS =-u0*sin(alpha0)*sin(lon);
  case 'const'
    uS=Param.uMax;
    vS=Param.vMax;
  case 'linear'
    uS=x(3);
    vS=0;
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
    