      # Model

      Grav=9.81;
cS=360;
Grav=9.81e0;
Cpd=1004.0e0;
Cvd=717.0e0;
Rd=Cpd-Cvd;
p0=1.0e5;
Cpv=1885.0e0;
Gamma=Cpd/Cvd;
kappa=Rd/Cpd;

# Sphere
Omega=2*pi/(24*3600);
RadEarth=6.37122e+6;

Deep=false;
HeightLimit=30000.0;
T0E=310.0;
T0P=240.0;
B=2.0;
K=3.0;
LapseRate=0.005;
U0=-0.5;
PertR=1.0/6.0;
Up=1.0;
PertExpR=0.1;
PertLon=pi/9.0;
PertLat=2.0*pi/9.0;
PertZ=15000.0;
      Lon=3.9269908169872414
      Lat=-0.6154797086703874
      Z=5250.0

  #     if Z>27500 && abs(Lat-pi/4)<=.1
  #       aa=3;
  #     end
      T0=0.5*(T0E+T0P);
      ConstA=1.0/LapseRate;
      ConstB=(T0-T0P)/(T0*T0P);
      ConstC=0.5*(K+2.0)*(T0E-T0P)/(T0E*T0P);
      ConstH=Rd*T0/Grav;
      ScaledZ=Z/(B*ConstH);
      Tau1=ConstA*LapseRate/T0*exp(LapseRate/T0*Z) +
        ConstB*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
      Tau2=ConstC*(1.0-2.0*ScaledZ*ScaledZ)*exp(-ScaledZ*ScaledZ);
      IntTau2=ConstC*Z*exp(-ScaledZ*ScaledZ);
      if Deep
        RRatio= R/EarthRadius;
      else
        RRatio = 1.0;
      end
      InteriorTerm=(RRatio*cos(Lat))^K -
        K/(K+2.0)*(RRatio*cos(Lat))^(K+2.0);
      Temperature=1.0/(RRatio*RRatio)/(Tau1-Tau2*InteriorTerm);

      InteriorTermU=(RRatio*cos(Lat))^(K-1.0) -
         (RRatio*cos(Lat))^(K+1.0);
      BigU=Grav/RadEarth*K *
        IntTau2*InteriorTermU*Temperature;
      if Deep
        RCosLat=R*cos(Lat);
      else
        RCosLat=RadEarth*cos(Lat);
      end
      OmegaRCosLat =Omega*RCosLat;

      # 		if (dOmegaRCosLat * dOmegaRCosLat + dRCosLat * dBigU < 0.0) {
      # 			_EXCEPTIONT("Negative discriminant detected.");
      # 		}

      @show Omega
      @show RCosLat
      uS=-OmegaRCosLat+sqrt(OmegaRCosLat*OmegaRCosLat+RCosLat*BigU);
      vS=0;
      # Exponential perturbation
      GreatCircleR=acos(sin(PertLat)*sin(Lat) +
        cos(PertLat)*cos(Lat)*cos(Lon-PertLon));

      GreatCircleR=GreatCircleR/PertExpR;

      # Tapered perturbation with height
      if Z < PertZ
        PertTaper=1.0-3.0*Z*Z/(PertZ * PertZ)+
           2.0*Z*Z*Z/(PertZ*PertZ*PertZ);
      else
        PertTaper=0.0;
      end

      # Apply perturbation in zonal velocity
      if GreatCircleR<1.0
        uSPert=Up*PertTaper*exp(-GreatCircleR*GreatCircleR);
      else
        uSPert=0.0;
      end
      uS=uS+uSPert;
