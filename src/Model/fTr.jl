function f = fTr(x,Param)
switch lower(Param.ProfTr)
  case 'bickley'
    f=Param.RhoTheta*sin(Param.k*x(2));
  case 'sphericalharmonics'
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    %f1=real(harmonicY(Param.lHar,Param.mHar,pi/2-lat,lon,r));
    f=harmonicS(Param.lHar,Param.mHar,lat,lon,r);
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


