function fVel(x,time,Global,Param)
    Model=Global.Model
    Phys=Global.Phys
    if Model.ProfVel == "SolidBody"
      uS=0;
      vS=0;
    elseif Model.ProfVel == "AdvectionSphereDCMIP"  
      (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
      Z=max(R-Phys.RadEarth,0);
      pZ = Phys.p0 * exp(-Z / Param.ScaleHeight)
      LonP = Lon - 2* pi * time / Param.tau
      k = 10 * Phys.RadEarth / Param.tau
      ua = k * sin(LonP)^2 * sin(2 * Lat) * cos(pi * time / Param.tau) +
       2 * pi * Phys.RadEarth / Param.tau * cos(Lat)
       ud = Param.omega_0 * Phys.RadEarth / Param.b / Param.p_top *
       cos(LonP) *
       cos(Lat)^2 *
       cos(2 * pi * time / Param.tau) *
       (-exp((pZ - Phys.p0) / Param.b / Param.p_top) + exp((Param.p_top - pZ) / Param.b / Param.p_top))
       uS = ua + ud
       vS = k * sin(2 * LonP) * cos(Lat) * cos(pi * time / Param.tau)
    elseif Model.ProfVel == "AdvectionTestDeform"  
      uS = -Param.uMax * (x[2] - Param.yC) * cospi(time / Param.EndTime)
      vS = Param.uMax * (x[1] - Param.xC) * cospi(time / Param.EndTime)
    elseif Model.ProfVel == "AdvectionSchaer"  
      vS = 0
      if x[3] <= Param.z1
        uS = 0  
      elseif x[3] <= Param.z2
        uS = Param.uMax * sin(0.5 * pi * (x[3] - Param.z1) / (Param.z2 - Param.z1))^2   
      else
        uS = Param.uMax      
      end  
      uS = Param.uMax      
    elseif Model.ProfVel == "SchaerSphere"
      (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
      uS = Param.uEq * cos(Lat)
      vS = 0.0
    elseif Model.ProfVel == "hyperdiffcart"
      if abs(x[1]-3000. * 1.e3)<=1000 * 1.e3 && abs(x[2])<=1000 * 1.e3
        uS=sqrt(5000);
        vS=sqrt(5000);
      else
        uS=0;
        vS=0;
      end
    elseif Model.ProfVel == "hyperdiff"
      (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3]);
      GreatCircleR=acos(sin(Param.PertLat)*sin(Lat)+
        cos(Param.PertLat)*cos(Lat)*cos(Lon-Param.PertLon));
      GreatCircleR=GreatCircleR/Param.PertExpR;
      # Apply perturbation in zonal velocity
      if GreatCircleR<6.0
        uS=sqrt(Param.ValDiff);
        vS=sqrt(Param.ValDiff);
      else
        uS=0.0;
        vS=0.0;
      end
    elseif Model.ProfVel == "BaroWaveSphere"
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
      IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
      if Param.Deep
        RRatio= R/Phys.EarthRadius;
      else
        RRatio = 1.0;
      end
      InteriorTerm=(RRatio*cos(Lat))^Param.K -
        Param.K/(Param.K+2.0)*(RRatio*cos(Lat))^(Param.K+2.0);
      Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);

      InteriorTermU=(RRatio*cos(Lat))^(Param.K-1.0) -
         (RRatio*cos(Lat))^(Param.K+1.0);
      BigU=Phys.Grav/Phys.RadEarth*Param.K *
        IntTau2*InteriorTermU*Temperature;
      if Param.Deep
        RCosLat=R*cos(Lat);
      else
        RCosLat=Phys.RadEarth*cos(Lat);
      end
      OmegaRCosLat =Phys.Omega*RCosLat;

      # 		if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
      # 			_EXCEPTIONT("Negative discriminant detected.");
      # 		}

      uS=-OmegaRCosLat+sqrt(OmegaRCosLat*OmegaRCosLat+RCosLat*BigU);
      vS=0;
      # Exponential perturbation
      GreatCircleR=acos(sin(Param.PertLat)*sin(Lat) +
        cos(Param.PertLat)*cos(Lat)*cos(Lon-Param.PertLon));

      GreatCircleR=GreatCircleR/Param.PertExpR;

      # Tapered perturbation with height
      if Z < Param.PertZ
        PertTaper=1.0-3.0*Z*Z/(Param.PertZ * Param.PertZ)+
           2.0*Z*Z*Z/(Param.PertZ*Param.PertZ*Param.PertZ);
      else
        PertTaper=0.0;
      end

      # Apply perturbation in zonal velocity
      if GreatCircleR<1.0
        uSPert=Param.Up*PertTaper*exp(-GreatCircleR*GreatCircleR);
      else
        uSPert=0.0;
      end
      uS=uS+uSPert;

    elseif Model.ProfVel == "barowavecart"
      eta=EtaFromZ(x[1],x[2],x[3],Param);
      uS=-Param.u0*sin(pi*x[2]/Param.Ly)^2*log(eta) *
        exp(-(log(eta)/Param.b)^2);
      if x[2]<0
        uS=-uS;
      else
        uS=uS+Param.uP*exp(-((x[1]-Param.xC)^2  +
          (x[2]-Param.yC)^2)/Param.Lp^2);
      end
      vS=0;
    elseif Model.ProfVel == "baldauf"
      uS=Param.uMax
      vS=0
    elseif Model.ProfVel == "barowavesphere"
  #     # parameters

    elseif Model.ProfVel == "gravityhill"
      uS=Param.uMax;
      vS=0;
    elseif Model.ProfVel == "sphericalharmonics"
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      (uS,vS)=harmonicVec(Param.lHar,Param.mHar,lat,lon,r);
    elseif Model.ProfVel == "galewsky"
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      if (lat<=Param.lat0G) || (lat>=Param.lat1G)
        uS=0;
      else
        uS=Param.uM/Param.eN*exp(1.0/((lat-Param.lat0G)*(lat-Param.lat1G)));
      end
      vS=0;
    elseif Model.ProfVel == "rossbyhaurwitz"
      omega6=7.8480e-6;
      K6=7.8480e-6;
      R6=4.0;
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      uS=r*omega6*cos(lat) + r*K6*(cos(lat))^(R6-1.0) *
        (R6*sin(lat)*sin(lat)-cos(lat)*cos(lat))*cos(R6*lon);
       vS=-r*K6*R6*(cos(lat))^(R6-1.0)*sin(lat)*sin(R6*lon);
    elseif Model.ProfVel == "hill"
      uMax=20;
      rotation_angle=0;
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      (rot_lon,rot_lat)=Rotate(lon,lat,rotation_angle);
      u_lat=uMax*cos(rot_lat);
      v_lat=0.0e0;
      if abs(rotation_angle)<1.0E-8
        uS = u_lat;
        vS = v_lat;
      else
        # rotate wind components
        (uS,vS) = turnwi(u_lat,v_lat,lon,lat,rot_lon,rot_lat,0.0e0,-0.5e0*pi+rotation_angle);
        if abs(uS)<1.0E-10
          uS=0.0e0;
        end
      end
    elseif Model.ProfVel == "bickley"
      U1 = cosh(x[2])^(-2);

      # ?? = exp(-(x2 + Param.l / 10)^2 / 2p.l^2) * cos(p.k * x1) * cos(p.k * x2)
      # Vortical velocity fields
      gaussian = exp(-(x[2] + Param.l / 10)^2 / (2*Param.l^2));
      u1 = gaussian * (x[2] + Param.l / 10) / Param.l^2 * cos(Param.k * x[1]) * cos(Param.k * x[2]);
      u1 = u1+Param.k * gaussian * cos(Param.k * x[1]) * sin(Param.k * x[2]);
      u2 = -Param.k * gaussian * sin(Param.k * x[1]) * cos(Param.k * x[2]);
      uS=U1+Param.Eps*u1;
      vS=Param.Eps*u2;
    elseif Model.ProfVel == "cosinebell"
      (lon,lat,r)=cart2sphere(x[1],x[2],x[3]);
      alpha0=pi/2;
      u0=2*pi*Phys.RadEarth/(86400*12);
      uS=u0*(cos(alpha0)*cos(lat)+sin(alpha0)*cos(lon)*sin(lat));
      vS =-u0*sin(alpha0)*sin(lon);
    elseif Model.ProfVel == "Const"
      uS=Param.uMax;
      vS=Param.vMax;
    elseif Model.ProfVel == "sin"
      uS=Param.uMax*sin(pi*x[1]/Param.Lx);
      vS=Param.vMax*sin(pi*x[2]/Param.Ly);
    elseif Model.ProfVel == "rand"
      pert = Param.pert * rand()
      uS=Param.uMax*pert;
      pert = Param.pert * rand()
      vS=Param.vMax*pert;
    elseif Model.ProfVel == "linear"
      uS=x[3];
      vS=0;
    else
      uS=0
      vS=0
    end
    return (uS,vS)
end

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

# zpih = pi*0.5d0;
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
# #-----------------------------------------------------------------------

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
