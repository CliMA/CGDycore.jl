function fRho(x,time,Phys,Global,Param,Profile)
  Model=Global.Model
  # global Omega uM lat0G lat1G eN
  str = lowercase(Model.ProfRho)
  if str == "solidbody"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    z=max(R-Phys.RadEarth,0);
    NBr=Param.NBr;
    Th0=Param.Th0;
    S=NBr*NBr/Phys.Grav;
    ThLoc=Th0*exp(z*S);
    pLoc=Phys.p0*(1-Phys.Grav/(Phys.Cpd*Param.Th0*S)*
        (1-exp(-S*z))).^(Phys.Cpd/Phys.Rd);
    Rho=pLoc./((pLoc/Phys.p0).^Phys.kappa*Phys.Rd.*ThLoc);
  elseif str == "testgrad"  
    Rho = x[1]^3+x[3]^3
  elseif str == "advectionspheredcmip"
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    Z=max(R-Phys.RadEarth,0);    
    p = Phys.p0 * exp(-Z / Param.ScaleHeight)
    Rho = p / Phys.Rd / Param.T_0
  elseif str == "advectionschaer"
    Rho = 1.0
  elseif str == "schaersphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Param.RadEarth,0);
    pLoc = Phys.p0 * exp(-Param.uEq * Param.uEq / (2.0 * Phys.Rd * Param.TEq) * sin(Lat)^2 - 
      Phys.Grav * Z / (Phys.Rd * Param.TEq))
    Rho = pLoc / (Phys.Rd * Param.TEq)
  elseif str == "hyperdiffcart"
    Rho=1;
  elseif str == "hyperdiff"
    Rho=1;
  elseif str == "heldsuarezcart"
    z = x[3]
    temp = Param.T_init + Param.lapse_rate * z 
    pres = Phys.p0 * (1 + Param.lapse_rate / Param.T_init * z)^(-Phys.Grav / Phys.Rd / Param.lapse_rate)
    Rho = pres / Phys.Rd / temp
  elseif str == "heldsuarezsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    z=max(R-Phys.RadEarth,0);
    temp = Param.T_init + Param.lapse_rate * z + rand() * 0.1 * (z < 5000)
    pres = Phys.p0 * (1 + Param.lapse_rate / Param.T_init * z)^(-Phys.Grav / Phys.Rd / Param.lapse_rate)
    Rho = pres / Phys.Rd / temp
  elseif str == "barowavecart"
    eta=EtaFromZ(x[1],x[2],x[3],Param);
    p=Phys.p0*eta;
    T=TBaroWave(x[1],x[2],eta,Param);
    Rho=p/(Phys.Rd*T);
elseif str == "barowavedrysphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Phys.Rd*T0/Phys.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z)+
    ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0)+
    ConstB*Z*exp(-ScaledZ*ScaledZ);
    IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
    if Param.Deep
      RRatio= R/Param.EarthRadius;
    else
      RRatio = 1.0;
    end
    InteriorTerm=(RRatio*cos(Lat))^Param.K -
    Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0);
    Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);
    Pressure=Phys.p0*exp(-Phys.Grav/Phys.Rd *
      (IntTau1-IntTau2*InteriorTerm));
    Rho=Pressure/(Phys.Rd*Temperature);
  elseif str == "baldaufsphere"
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
    r=r-Phys.RadEarth / Param.ScaleRad
    p=Phys.p0*exp(-Phys.Grav*r/(Phys.Rd*Param.T0));
    Rho=p/(Phys.Rd*Param.T0);
    T=p/(Rho*Phys.Rd);
    ExpDist = Phys.Grav / (Phys.Rd * Param.T0)
    T=T+Param.DeltaT*exp(-ExpDist*r)*exp(Param.eta * (sin(lat) -1.0))*sin(pi*r/Param.H);
    Rho=p/(Phys.Rd*T);
  elseif str == "barowavecart"
    eta=EtaFromZ(x[1],x[2],x[3],Param);
    p=Phys.p0*eta;
    T=TBaroWave(x[1],x[2],eta,Param);
    Rho=p/(Phys.Rd*T);
  elseif str == "barowavemoistsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Phys.Rd*T0/Phys.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z)+
    ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0)+
      ConstB*Z*exp(-ScaledZ*ScaledZ);
    IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
    if Param.Deep
      RRatio= R/Phys.EarthRadius;
    else
      RRatio = 1.0;
    end
    InteriorTerm=(RRatio*cos(Lat))^Param.K -
      Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0);
    Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);
    Pressure=Phys.p0*exp(-Phys.Grav/Phys.Rd *
        (IntTau1-IntTau2*InteriorTerm));
    Rho=Pressure/(Phys.Rd*Temperature)
  elseif str == "baldaufsphere"
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      r=r-Phys.RadEarth;
      p=Phys.p0*exp(-Phys.Grav*r/(Phys.Rd*Param.T0));
      Rho=p/(Phys.Rd*Param.T0);
      T=p/(Rho*Phys.Rd);
      d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
      T=T+Param.DeltaT*exp(-Param.ExpDist*d)*sin(pi*r/Param.H);
      Rho=p/(Phys.Rd*T);
  elseif str == "baldaufcart"
      delta=Phys.Grav/(Phys.Rd*Param.T0);

      p=Phys.p0*exp(-delta*x[3]);

      dT=Param.DeltaT*exp(-(x[2]-Param.yc)^2/Param.d^2)*sin(pi*x[3]/Param.H);

      TLoc=Param.T0+exp(delta/2*x[3])*dT;
      Rho=p/(Phys.Rd*TLoc);
  elseif str == "bryanfritsch"
    z = x[3]
    @views zP = Profile[:,1]
    iz = 1000
    for i = 2:size(zP,1)
      if z <= zP[i]
        iz = i - 1
        break
      end  
    end
    z_l = zP[iz]
    Rho_l = Profile[iz,2]
    Theta_l = Profile[iz,3]
    RhoV_l = Profile[iz,4]
    RhoC_l = Profile[iz,5]
    z_r = zP[iz+1]
    Rho_r = Profile[iz+1,2]
    Theta_r = Profile[iz+1,3]
    RhoV_r = Profile[iz+1,4]
    RhoC_r = Profile[iz+1,5]
    Rho = (Rho_r * (z - z_l) + Rho_l * (z_r - z)) / (z_r - z_l)
    Theta = (Theta_r * (z - z_l) + Theta_l * (z_r - z)) / (z_r - z_l)
    RhoV = (RhoV_r * (z - z_l) + RhoV_l * (z_r - z)) / (z_r - z_l)
    RhoC = (RhoC_r * (z - z_l) + RhoC_l * (z_r - z)) / (z_r - z_l)

    Rho, Theta, RhoV, RhoC = PerturbMoistProfile(x, Rho, Rho*Theta, RhoV, RhoC, Phys, Param)

  elseif str == "warmbubble2dx"
      Grav=Phys.Grav;
      p0=Phys.p0;
      Rd=Phys.Rd;
      kappa=Phys.kappa;
      Th0=Param.Th0;
      DeltaTh=Param.DeltaTh;
      xC0=Param.xC0;
      zC0=Param.zC0;
      rC0=Param.rC0;
      x3=x[3];
      x1=x[1];
      pLoc=p0*(1-Grav*x3*kappa/(Rd*Th0))^(1/kappa);
      rr=sqrt((x1-xC0)^2+(x3-zC0)^2);
      ThLoc=Th0;
      if rr<rC0
        ThLoc=ThLoc+DeltaTh*cos(0.5*pi*rr/rC0)^2;
      end
      Rho=pLoc/((pLoc/p0)^kappa*Rd.*ThLoc);
  elseif str == "warmbubble2dy"
      Grav=Phys.Grav;
      p0=Phys.p0;
      Rd=Phys.Rd;
      kappa=Phys.kappa;
      Th0=Param.Th0;
      DeltaTh=Param.DeltaTh;
      yC0=Param.yC0;
      zC0=Param.zC0;
      rC0=Param.rC0;
      x3=x[3];
      x2=x[2];
      pLoc=p0*(1-Grav*x3*kappa/(Rd*Th0))^(1/kappa);
      rr=sqrt((x2-yC0)^2+(x3-zC0)^2);
      ThLoc=Th0;
      if rr<rC0
        ThLoc=ThLoc+DeltaTh*cos(0.5*pi*rr/rC0)^2;
      end
      Rho=pLoc/((pLoc/p0)^kappa*Rd.*ThLoc);
  elseif str == "densitycurrent"
      T0=Param.T0;
      DeltaT=Param.DeltaT;
      xC0=Param.xC0;
      zC0=Param.zC0;
      xrC0=Param.xrC0;
      zrC0=Param.zrC0;
      x3=x[3];
      x1=x[1];
      pLoc=Phys.p0*(1-Phys.kappa*Phys.Grav*x[3] /
        (Phys.Rd*T0))^(Phys.Cpd/Phys.Rd);
      Rad=sqrt(((x1-xC0)/xrC0)^2+((x3-zC0)/zrC0)^2);
      ThLoc=T0;
      if Rad<1.0e0
        ThLoc=ThLoc+DeltaT*(cos(pi*Rad)+1.0)/2.0*(pLoc/Phys.p0)^(-Phys.kappa);
      end
      Rho=pLoc/((pLoc/Phys.p0)^Phys.kappa*Phys.Rd*ThLoc);
  elseif str == "gravityhill" || str == "schaercart" || str == "agnesicart"
      z=x[3];
      NBr=Param.NBr;
      Grav=Phys.Grav;
      p0=Phys.p0;
      Cpd=Phys.Cpd;
      Rd=Phys.Rd;
      kappa=Phys.kappa;
      Th0=Param.Th0;
      S=NBr*NBr/Grav;
      ThLoc=Th0*exp(z*S);
      pLoc=p0*(1.0-Grav/(Cpd*Th0*S)*(1.0-exp(-S*z))).^(Cpd/Rd);
      Rho=pLoc/((pLoc/p0)^kappa*Rd*ThLoc);
  elseif str == "inertiagravitywave"
      z=x[3];
      NBr=Param.NBr;
      Grav=Phys.Grav;
      Th0=Param.Th0;
      S=NBr*NBr/Grav;
      ThB=Th0*exp(z*S);
      pLoc=Phys.p0*(1-Phys.Grav/(Phys.Cpd*Th0*S)*(1-exp(-S*z))).^(Phys.Cpd/Phys.Rd);
      ThLoc=ThB+Param.DeltaTh*sin(pi*x[3]/Param.H)./(1+(x[2]-Param.yC).^2/Param.a^2);
      Rho=pLoc./((pLoc/Phys.p0).^Phys.kappa*Phys.Rd.*ThLoc);
  elseif str == "galewski"
      Grav=Phys.Grav;
      Omega=Phys.Omega;
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      r=Phys.RadEarth;
      Rho=(Grav*Param.H0G-(simpson(-0.5*pi,lat,r,pi/100.0,integrandG,Param)))/Grav +
        Param.hH*cos(lat)*exp(-((lon-pi)/Param.alphaG)^2.0)*exp(-((pi/4.0-lat)/Param.betaG)^2.0);
  elseif str == "rossbyhaurwitz"
      Grav=Phys.Grav;
      Omega=Phys.Omega;
      H06=8000.0;
      omega6=7.8480e-6;
      K6=7.8480e-6;
      R6=4.0;
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      A=0.5*omega6*(2.0*Omega+omega6)*cos(lat)*cos(lat) +
        0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat) +
        (2.0*R6*R6-R6-2.0)-2.0*R6*R6*(cos(lat))^(-2));
      B=2.0*(Omega+omega6)*K6/(R6+1.0)/(R6+2.0)*
        (cos(lat))^R6*((R6*R6+2.0*R6+2.0)-
        (R6+1.0)*(R6+1.0)*cos(lat)*cos(lat));
      C=0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat)-(R6+2.0));
      Rho=(Grav*H06+r*r*(A+B*cos(R6*lon)+C*cos(2.0*R6*lon)))/Grav;
  elseif str == "hill"
      lon0=3.0e0/2.0e0*pi;
      lat0=pi/6.0e0;
      H05=5960.0e0;
      hS=2000.0e0;
      RadiusC=pi/9.0e0;
      UMax=20.0e0;
      rotation_angle=0;
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      (rot_lon,rot_lat)=Rotate(lon,lat,rotation_angle);
      (rot_lon0,rot_lat0)=Rotate(lon0,lat0,rotation_angle);
      r=sqrt(min(RadiusC*RadiusC,
        (rot_lon-rot_lon0)*(rot_lon-rot_lon0)+(rot_lat-rot_lat0)*(rot_lat-rot_lat0)));
      HeightLoc=hS*(1-r/RadiusC);
      Rho=(Phys.Grav*H05-(Phys.RadEarth*Phys.Omega*UMax+0.5*UMax*UMax)*sin(rot_lat)*sin(rot_lat))/
        Phys.Grav-HeightLoc;
  elseif str == "spherical"
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
      if abs(d)<=0.8
        Rho=100*cos(pi*d/0.8/2)^2+100;
      else
        Rho=100.0;
      end
  elseif str == "cosinebell"
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      r1=Phys.RadEarth*acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
      R3=Phys.RadEarth/3.0e0;
      if r1<=R3
        Rho=1000.0e0/2.0e0*(1.0e0+cos(pi*r1/R3));
      else
        Rho=0.0e0;
      end
  elseif str == "sphericalsmooth"
      Rho=exp(-((-x*normal-1.0)/0.1)^2);
      #Rho=exp(-((-x(2)-1.0)/0.4)^2)
  elseif str == "quad"
      if abs(x[1]-xM[1])<xH[1] && abs(x[2]-xM[2])<xH[2]
        Rho=1;
      else
        Rho=0;
      end
  elseif str == "cosine"
      d=sqrt((x[1]-xM[1])^2+(x[2]-xM[2])^2);
      if d<=rH
        Rho=cos(pi*d/rH/2)^2;
      else
        Rho=0.0;
      end
  elseif str == "constant"
      Rho=1;
  elseif str == "linear"
      Rho=x[1]+1;
  elseif str == "bickley"
      Rho=Param.RhoTheta;
  elseif str == "isothermalcart"
    pLoc = Phys.p0 * exp(-Phys.Grav * x[3] / (Phys.Rd * Param.TEq))
    Rho = pLoc / (Phys.Rd * Param.TEq)
  elseif str == "isothermalsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    pLoc = Phys.p0 * exp(-Phys.Grav * Z / (Phys.Rd * Param.TEq))
    Rho = pLoc / (Phys.Rd * Param.TEq)  
  elseif str == "decayingtemperatureprofile"
    H_sfc = Phys.Rd * Param.T_virt_surf / Phys.Grav
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
    z=r-Phys.RadEarth;
    zprime = z / Param.H_t
    tanh_zprime = tanh(zprime)
    Delta_Tv = Param.T_virt_surf - Param.T_min_ref
    Tv = Param.T_virt_surf - Delta_Tv * tanh_zprime
    Delta_Tvprime = Delta_Tv / Param.T_virt_surf
    p = -Param.H_t * (zprime + Delta_Tvprime * (log(1 - Delta_Tvprime * tanh_zprime) - log(1 + tanh_zprime) + zprime))
    p /= H_sfc * (1 - Delta_Tvprime^2)
    p = Phys.p0 * exp(p)
    Rho = p / (Phys.Rd * Tv)
  else
    Rho=1.0  
  end
  return Rho
end

function integrandG(tau,RadEarth,Param)
  f=2.0*Param.Omega*sin(tau);
  if (tau<=Param.lat0G) || (tau>=Param.lat1G)
    uStart=0.0;
  else
    uStart=Param.uM/Param.eN*exp(1.0/((tau-Param.lat0G)*(tau-Param.lat1G)));
  end
  if abs(tau)<0.5*pi
    intG=(RadEarth*f+uStart*tan(tau))*uStart;
  else
    intG=0.0;
  end
  return intG
end

function Rotate(lon,lat,rotation_angle)
  if abs(rotation_angle)<1.0e-8
    rot_lon = lon;
    rot_lat = lat;
  else
    (rot_lon,rot_lat)=regrot(lon,lat,0.0e0,-0.5e0*pi+rotation_angle);
  end
  return rot_lon,rot_lat
end

function regrot(pxreg,pyreg,pxcen,pycen)

#----------------------------------------------------------------------
#
#*  conversion between regular and rotated spherical coordinates.
#*
#*  pxreg     longitudes of the regular coordinates
#*  pyreg     latitudes of the regular coordinates
#*  pxrot     longitudes of the rotated coordinates
#*  pyrot     latitudes of the rotated coordinates
#*            all coordinates given in degrees n (negative for s)
#*            and degrees e (negative values for w)
#*  pxcen     regular longitude of the south pole of the rotated grid
#*  pycen     regular latitude of the south pole of the rotated grid
#*
#*  kcall=-1: find regular as functions of rotated coordinates.
#*  kcall= 1: find rotated as functions of regular coordinates.
#
#-----------------------------------------------------------------------
#

zpih = pi*0.5e0;
zsycen = sin((pycen+zpih));
zcycen = cos((pycen+zpih));

zxmxc= pxreg - pxcen;
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
return (pxrot,pyrot)
end

function turnwi(puarg,pvarg,
                pxreg,pyreg,pxrot,pyrot,
                pxcen,pycen)
#
#-----------------------------------------------------------------------
#
#*  turn horizontal velocity components between regular and
#*  rotated spherical coordinates.
#
#*  puarg : input u components
#*  pvarg : input v components
#*  pures : output u components
#*  pvres : output v components
#*  pa    : transformation coefficients
#*  pb    :    -"-
#*  pc    :    -"-
#*  pd    :    -"-
#*  pxreg : regular longitudes
#*  pyreg : regular latitudes
#*  pxrot : rotated longitudes
#*  pyrot : rotated latitudes
#*  kxdim              : dimension in the x (longitude) direction
#*  kydim              : dimension in the y (latitude) direction
#*  kx                 : number of gridpoints in the x direction
#*  ky                 : number of gridpoints in the y direction
#*  pxcen              : regular longitude of the south pole of the
#*                       transformed grid
#*  pycen              : regular latitude of the south pole of the
#*                       transformed grid
#*
#*  kcall < 0          : find wind components in regular coordinates
#*                       from wind components in rotated coordinates
#*  kcall > 0          : find wind components in rotated coordinates
#*                       from wind components in regular coordinates
#*  note that all coordinates are given in degrees n and degrees e.
#*     (negative values for s and w)
#
#-----------------------------------------------------------------------

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
   pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -
        zcxmxc*zsxrot*zsyrot;
   pc =-zsyc*zsxrot/zcyreg;
   pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg;

   pures = pa*puarg + pb*pvarg;
   pvres = pc*puarg + pd*pvarg;
   return pures,pvres
end

