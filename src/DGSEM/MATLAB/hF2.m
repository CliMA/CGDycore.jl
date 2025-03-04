function h = hF2(x,Param)
global lat0 lon0 xM xH rH normal uM lat0G lat1G eN Omega
switch lower(Param.Profh)
  case 'galewsky'
    Grav=Param.Grav;
    Omega=Param.Omega;
    alphaG=1.0/3.0;
    betaG=1.0/15.0;
    hH=120.0;
    H0G=10000.0;
    uM=80.0;
    lat0G=pi/7.0;
    lat1G=pi/2.0-lat0G;
    eN=exp(-4.0/(lat1G-lat0G)^2.0);
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    h=(Grav*H0G-(simpson(-0.5*pi,lat,r,pi/100.0,@integrandG)))/Grav...
    +hH*cos(lat)*exp(-((lon-pi)/alphaG)^2.0)*exp(-((pi/4.0-lat)/betaG)^2.0);
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
    h=(Grav*H06+r*r*(A+B*cos(R6*lon)+C*cos(2.0*R6*lon)))/Grav;
  case 'spherical'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    d=acos(sin(lat0)*sin(lat)+cos(lat0)*cos(lat)*cos(lon-lon0));
    if abs(d)<=0.8
      h=100*cos(pi*d/0.8/2)^2+100;
    else
      h=100.0;
    end
  case 'sphericalsmooth'
    h=exp(-((-x*normal-1.0)/0.1)^2);
    %h=exp(-((-x(2)-1.0)/0.4)^2)
  case 'quad'
    if abs(x(1)-xM(1))<xH(1) && abs(x(2)-xM(2))<xH(2)
      h=1;
    else
      h=0;
    end
  case 'cosine'
    d=sqrt((x(1)-xM(1))^2+(x(2)-xM(2))^2);
    if d<=rH
      h=cos(pi*d/rH/2)^2;
    else
      h=0.0;
    end
  case 'constant'
    h=1;
  case 'linear'
    h=x(1)+1;
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

function res=simpson(x0,xN,r,dx,f)
  n=(xN-x0)/dx+1;
  h=(xN-x0)/n;
  res=0.5*(f(x0,r)+f(xN,r));
  xi=x0;
  for i=1:(n-1)
    xi=xi+h;
    res=res+f(xi,r);
  end
  xi=x0-0.5*h;
  for i=1:n
    xi=xi+h;
    res=res+2.0*f(xi,r);
  end
  res=res*h/3.0;
end


