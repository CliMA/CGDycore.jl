function Ros=MethodRos(name)
switch name
  case 'RosRK3'
    Ros.nStage=3;
    Ros.beta0=1.0e0;
    beta0=Ros.beta0;
    %     beta0=1.0e0/2.0e0
    %    beta0=(3.0e0+SQRT(3.0e0))/6.0e0
    Ros.alpha(1,1)=1.0e0/(3.0e0*beta0);
    Ros.alpha(2,1)=(-1.0e0 + 12.0e0*beta0^2.0e0)/(18.0e0*beta0^2.0e0*(-1.0e0+4.0e0*beta0));
    Ros.alpha(2,2)=1.0e0/(2.0e0*beta0);
    Ros.alpha(3,1)=(1.0e0-3.0e0*beta0*(7.0e0+16.0e0*beta0*(-2.0e0+3.0e0*beta0)))...
      /(36.0e0*beta0^3.0e0*(-1.0e0+4.0e0*beta0));
    Ros.alpha(3,2)=(-1.0e0+12.0e0*beta0)/(4.0e0*beta0^2.0e0);
    Ros.alpha(3,3)=1.0e0/beta0;
    Ros.beta(2,1)=(1.0e0-12.0e0*beta0^2.0e0)/(beta0^2.0e0*(-9.0e0+36.0e0*beta0));
    Ros.beta(3,1)=(-1.0e0+3.0e0*beta0*(7.0e0+16.0e0*beta0*(-2.0e0+3.0e0*beta0)))...
      /(36.*beta0^3.0e0*(-1.0e0+4.0e0*beta0));
    Ros.beta(3,2)=(0.250e0-3.0e0*beta0)/beta0^2.0e0;
    
    Ros.alpha=Ros.alpha*Ros.beta0;
    Ros.beta=Ros.beta*Ros.beta0; 
  case 'RosEul'
    Ros.nStage=1;
    Ros.beta0=1.0e0;
    Ros.alpha(1,1)=1.0e0;
end
end
