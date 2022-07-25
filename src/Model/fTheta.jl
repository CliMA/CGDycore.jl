function fTheta(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.ProfTheta)
  if str == "solidbody"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    z=max(R-Phys.RadEarth,0);
    NBr=Param.NBr;
    Grav=Phys.Grav;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    Th=Th0*exp(z*S);
  elseif str == "schaersphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Param.RadEarth,0);
    pLoc = Phys.p0 * exp(-Param.uEq * Param.uEq / (2.0 * Phys.Rd * Param.TEq) * sin(Lat)^2 -
      Phys.Grav * Z / (Phys.Rd * Param.TEq))
    Th=Param.TEq * (Phys.p0 / pLoc)^(Phys.Rd / Phys.Cpd)
  elseif str == "hyperdiffcart"
    if abs(x[1]-3000. * 1.e3)<=1000 * 1.e3 && abs(x[2])<=1000 * 1.e3
      Th=100.;
    else
      Th=0;
    end
  elseif str == "hyperdiff"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    GreatCircleR=acos(sin(Param.PertLat)*sin(Lat) +
      cos(Param.PertLat)*cos(Lat)*cos(Lon-Param.PertLon));
    GreatCircleR=GreatCircleR/Param.PertExpR;
    # Apply perturbation in zonal velocity
    if GreatCircleR<6.0
      Th=sqrt(2*Param.ValDiff);
    else
      Th=0.0;
    end
  elseif str == "heldsuarezcart"
        z = x[3]
        pert = Param.pert*rand() 
        temp = Param.T_init + Param.lapse_rate * z + pert * 0.1 * (z < 5000)
        pres = Phys.p0 * (1 + Param.lapse_rate / Param.T_init * z)^(-Phys.Grav / Phys.Rd / Param.lapse_rate)
        Th = temp * (Phys.p0 / pres)^Phys.kappa
  elseif str == "heldsuarezsphere"
        (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
        z=max(R-Phys.RadEarth,0);
        temp = Param.T_init + Param.lapse_rate * z + rand() * 0.1 * (z < 5000)
        pres = Phys.p0 * (1 + Param.lapse_rate / Param.T_init * z)^(-Phys.Grav / Phys.Rd / Param.lapse_rate)
        Th = temp * (Phys.p0 / pres)^Phys.kappa
  elseif str == "barowavecart"
    eta=EtaFromZ(x[1],x[2],x[3],Param);
    p=Phys.p0*eta;
    T=TBaroWave(x[1],x[2],eta,Param);
    Th=T*(Phys.p0/p)^(Phys.Rd/Phys.Cpd);
  elseif str == "barowavesphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Phys.Rd*T0/Phys.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z) +
      ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0) +
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
    Th=Temperature*(Phys.p0/Pressure)^(Phys.Rd/Phys.Cpd);
  elseif str == "barowavemoistsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
    Z=max(R-Phys.RadEarth,0);
    T0=0.5*(Param.T0E+Param.T0P);
    ConstA=1.0/Param.LapseRate;
    ConstB=(T0-Param.T0P)/(T0*Param.T0P);
    ConstC=0.5*(Param.K+2.0)*(Param.T0E-Param.T0P)/(Param.T0E*Param.T0P);
    ConstH=Phys.Rd*T0/Phys.Grav;
    ScaledZ=Z/(Param.B*ConstH);
    Tau1=ConstA*Param.LapseRate/T0*exp(Param.LapseRate/T0*Z) +
      ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
    IntTau1=ConstA*(exp(Param.LapseRate/T0*Z)-1.0) +
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
    Rho=Pressure/(Phys.Rd*Temperature)
    q_0 = Param.q_0                # Maximum specific humidity (default: 0.018)
    q_t = Param.q_t                # Specific humidity above artificial tropopause
    lat_w = 2.0*pi / 9.0
    p_w = 34.0e3
    eta_crit = p_w / Phys.p0
    eta = Pressure/Phys.p0
    if eta>eta_crit
      qv = q_0*exp(-(Lat/lat_w)^4)*exp(-((eta-1.)*Phys.p0/p_w)^2)
    else
      qv = q_t
    end
    Mvap = 0.608e0               #  Ratio of molar mass of dry air/water
    Temperature = Temperature / (1.0 + Mvap * qv)

    #!-----------------------------------------------------
    #!   Initialize virtual potential temperature
    #!-----------------------------------------------------
    #thetav = t * (1.d0 + 0.61d0 * q) * (p0 / p)**(Rd / cp)

    RhoV=Rho * qv
    RhoD = Rho - RhoV 
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV
    Rm  = Phys.Rd * RhoD + Phys.Rv * RhoV
    Th=Temperature*(Phys.p0/Pressure)^(Rm/Cpml) * (RhoD + (Phys.Rv / Phys.Rd) *RhoV) / (RhoD + RhoV)
  elseif str == "baldaufsphere"
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
    r=r-Phys.RadEarth;
    p=Phys.p0*exp(-Phys.Grav*r/(Phys.Rd*Param.T0));
    Rho=p/(Phys.Rd*Param.T0);
    T=p/(Rho*Phys.Rd);
    d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    T=T+Param.DeltaT*exp(-Param.ExpDist*d)*sin(pi*r/Param.H);
    Th=T*(Phys.p0/p)^(Phys.Rd/Phys.Cpd);
  elseif str == "baldaufcart"
    delta=Phys.Grav/(Phys.Rd*Param.T0);
    p=Phys.p0*exp(-delta*x[3]);

    dT=Param.DeltaT*exp(-(x[1]-Param.xc)^2/Param.d^2)*sin(pi*x[3]/Param.H);

    TLoc=Param.T0+exp(delta/2*x[3])*dT;
    Th=TLoc*(Phys.p0/p)^(Phys.Rd/Phys.Cpd);
  elseif str == "warmbubble2d"
    Th0=Param.Th0;
    DeltaTh=Param.DeltaTh;
    xC0=Param.xC0;
    zC0=Param.zC0;
    rC0=Param.rC0;
    x3=x[3];
    x1=x[1];
    rr=sqrt((x1-xC0)^2+(x3-zC0)^2);
    ThLoc=Th0;
    if rr<rC0
      ThLoc=ThLoc+DeltaTh*cos(0.5*pi*rr/rC0)^2;
    end
    Th=ThLoc;
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
    Th=T0;
    if Rad<1.0e0
      Th=Th+DeltaT*(cos(pi*Rad)+1.0)/2.0*(pLoc/Phys.p0)^(-Phys.kappa);
    end
  elseif str == "gravityhill" || str == "schaercart"
    z=x[3];
    NBr=Param.NBr;
    Grav=Phys.Grav;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    Th=Th0*exp(z*S);
  elseif str == "inertiagravitywave"
    z=x[3];
    NBr=Param.NBr;
    Grav=Phys.Grav;
    Th0=Param.Th0;
    S=NBr*NBr/Grav;
    ThB=Th0*exp(z*S);
    Th=ThB+Param.DeltaTh*sin(pi*x[3]/Param.H)/(1+(x[2]-Param.yC)^2/Param.a^2);
  elseif str == "galewsky"
    Th=1;
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
    B=2.0*(Omega+omega6)*K6/(R6+1.0)/(R6+2.0) *
      (cos(lat))^R6*((R6*R6+2.0*R6+2.0) -
        (R6+1.0)*(R6+1.0)*cos(lat)*cos(lat));
    C=0.25*K6*K6*(cos(lat))^(2.0*R6)*((R6+1.0)*cos(lat)*cos(lat)-(R6+2.0));
    Th=(Grav*H06+r*r*(A+B*cos(R6*lon)+C*cos(2.0*R6*lon)))/Grav;
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
    Th=(Phys.Grav*H05-(Phys.RadEarth*Phys.Omega*UMax+0.5*UMax*UMax)*sin(rot_lat)*sin(rot_lat)) /Phys.Grav-HeightLoc;
  elseif str == "spherical"
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
    d=acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    if abs(d)<=0.8
      Th=100*cos(pi*d/0.8/2)^2+100;
    else
      Th=100.0;
    end
  elseif str == "cosinebell"
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
    r1=Phys.RadEarth*acos(sin(Param.lat0)*sin(lat)+cos(Param.lat0)*cos(lat)*cos(lon-Param.lon0));
    R3=Phys.RadEarth/3.0e0;
    if r1<=R3
      Th=1000.0e0/2.0e0*(1.0e0+cos(pi*r1/R3));
    else
      Th=0.0e0;
    end
  elseif str == "sphericalsmooth"
    Th=exp(-((-x*normal-1.0)/0.1)^2);
    #Th=exp(-((-x[2]-1.0)/0.4)^2)
  elseif str == "quad"
    if abs(x[1]-xM[1])<xH[1] && abs(x[2]-xM[2])<xH[2]
      Th=1;
    else
      Th=0;
    end
  elseif str == "cosine"
    d=sqrt((x[1]-xM[1])^2+(x[2]-xM[2])^2);
    if d<=rH
      Th=cos(pi*d/rH/2)^2;
    else
      Th=0.0;
    end
  elseif str == "constant"
    Th=1;
  elseif str == "linear"
    Th=x[1]+1;
  elseif str == "bickley"
    Th=Param.RhoTheta;
  elseif str == "isothermal"
    pLoc = Phys.p0 * exp(-Phys.Grav * x[3] / (Phys.Rd * Param.TEq))
    Th=Param.TEq * (Phys.p0 / pLoc)^(Phys.Rd / Phys.Cpd)  
  elseif str == "decayingtemperatureprofile"  
    H_sfc = Phys.Rd * Param.T_virt_surf / Phys.Grav
    (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
    z=r-Phys.RadEarth;
    zprime = z / Param.H_t 
    tanh_zprime = tanh(zprime)
    Delta_Tv = Param.T_virt_surf - Param.T_min_ref
    Tv = Param.T_virt_surf - Delta_Tv * tanh_zprime + rand() * 0.1 * (z < 5000)
    Delta_Tvprime = Delta_Tv / Param.T_virt_surf
    p = -Param.H_t * (zprime + Delta_Tvprime * (log(1 - Delta_Tvprime * tanh_zprime) - log(1 + tanh_zprime) + zprime))
    p /= H_sfc * (1 - Delta_Tvprime^2)
    p = Phys.p0 * exp(p)
    Th=Tv * (Phys.p0 / p)^(Phys.Rd / Phys.Cpd)  
  end
  return Th
end

# function intG=integrandG(tau,RadEarth)
# global Omega uM lat0G lat1G eN
# f=2.0*Omega*sin(tau);
# if (tau<=lat0G) || (tau>=lat1G)
#   uStart=0.0;
# else
#   uStart=uM/eN*exp(1.0/((tau-lat0G)*(tau-lat1G)));
# end
# if abs(tau)<0.5*pi
#   intG=(RadEarth*f+uStart*tan(tau))*uStart;
# else
#   intG=0.0;
# end
# end

# function [rot_lon,rot_lat]=Rotate(lon,lat,rotation_angle)
# if abs(rotation_angle)<1.0e-8
#   rot_lon = lon;
#   rot_lat = lat;
# else
#   [rot_lon,rot_lat]=regrot(lon,lat,0.0e0,-0.5e0*pi+rotation_angle);
# end
# end

# function [pxrot,pyrot]=regrot(pxreg,pyreg,pxcen,pycen)

# #----------------------------------------------------------------------
# #
# #*    conversion between regular and rotated spherical coordinates.
# #*
# #*    pxreg     longitudes of the regular coordinates
# #*    pyreg     latitudes of the regular coordinates
# #*    pxrot     longitudes of the rotated coordinates
# #*    pyrot     latitudes of the rotated coordinates
# #*              all coordinates given in degrees n (negative for s)
# #*              and degrees e (negative values for w)
# #*    pxcen     regular longitude of the south pole of the rotated grid
# #*    pycen     regular latitude of the south pole of the rotated grid
# #*
# #*    kcall=-1: find regular as functions of rotated coordinates.
# #*    kcall= 1: find rotated as functions of regular coordinates.
# #
# #-----------------------------------------------------------------------
# #

# zpih = pi*0.5e0;
# zsycen = sin((pycen+zpih));
# zcycen = cos((pycen+zpih));

# zxmxc  = pxreg - pxcen;
# zsxmxc = sin(zxmxc);
# zcxmxc = cos(zxmxc);
# zsyreg = sin(pyreg);
# zcyreg = cos(pyreg);
# zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc;
# zsyrot = max(zsyrot,-1.e0);
# zsyrot = min(zsyrot,+1.e0);

# pyrot = asin(zsyrot);

# zcyrot = cos(pyrot);
# zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot;
# zcxrot = max(zcxrot,-1.e0);
# zcxrot = min(zcxrot,+1.e0);
# zsxrot = zcyreg*zsxmxc/zcyrot;

# pxrot = acos(zcxrot);

# if zsxrot<0.0
#   pxrot = -pxrot;
# end
# end

# function [pures,pvres]=turnwi(puarg,pvarg, 
#                   pxreg,pyreg,pxrot,pyrot, 
#                   pxcen,pycen)
# #
# #-----------------------------------------------------------------------
# #
# #*    turn horizontal velocity components between regular and
# #*    rotated spherical coordinates.
# #
# #*    puarg : input u components
# #*    pvarg : input v components
# #*    pures : output u components
# #*    pvres : output v components
# #*    pa    : transformation coefficients
# #*    pb    :    -"-
# #*    pc    :    -"-
# #*    pd    :    -"-
# #*    pxreg : regular longitudes
# #*    pyreg : regular latitudes
# #*    pxrot : rotated longitudes
# #*    pyrot : rotated latitudes
# #*    kxdim              : dimension in the x (longitude) direction
# #*    kydim              : dimension in the y (latitude) direction
# #*    kx                 : number of gridpoints in the x direction
# #*    ky                 : number of gridpoints in the y direction
# #*    pxcen              : regular longitude of the south pole of the
# #*                         transformed grid
# #*    pycen              : regular latitude of the south pole of the
# #*                         transformed grid
# #*
# #*    kcall < 0          : find wind components in regular coordinates
# #*                         from wind components in rotated coordinates
# #*    kcall > 0          : find wind components in rotated coordinates
# #*                         from wind components in regular coordinates
# #*    note that all coordinates are given in degrees n and degrees e.
# #*       (negative values for s and w)
# #
# !-----------------------------------------------------------------------

#      zpih = pi*0.5e0;
#      zsyc = sin(pycen+zpih);
#      zcyc = cos(pycen+zpih);

#      zsxreg = sin(pxreg);
#      zcxreg = cos(pxreg);
#      zsyreg = sin(pyreg);
#      zcyreg = cos(pyreg);

#      zxmxc  = pxreg - pxcen;
#      zsxmxc = sin(zxmxc);
#      zcxmxc = cos(zxmxc);

#      zsxrot = sin(pxrot);
#      zcxrot = cos(pxrot);
#      zsyrot = sin(pyrot);
#      zcyrot = cos(pyrot);

#      pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot;
#      pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot - 
#           zcxmxc*zsxrot*zsyrot;
#      pc =-zsyc*zsxrot/zcyreg;
#      pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg;
#      !
#      pures = pa*puarg + pb*pvarg;
#      pvres = pc*puarg + pd*pvarg;
# end
