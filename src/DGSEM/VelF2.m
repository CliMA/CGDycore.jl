function VelSp = VelF2(x,Param)
switch lower(Param.ProfuVel)
  case 'galewsky'
    uM=80.0;
    lat0G=pi/7.0;
    lat1G=pi/2.0-lat0G;
    eN=exp(-4.0/(lat1G-lat0G)^2.0);
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    VelSp=zeros(2,1);
    if (lat<=lat0G) || (lat>=lat1G)
      VelSp(1)=0;
    else
      VelSp(1)=uM/eN*exp(1.0/((lat-lat0G)*(lat-lat1G)));
    end
    VelSp(2)=0;
  case 'rossbyhaurwitz'
    omega6=7.8480e-6; 
    K6=7.8480e-6;
    R6=4.0;
    [lon,lat,r]=cart2sphere(x(1),x(2),x(3));
    VelSp=zeros(2,1);
    VelSp(1)=r*omega6*cos(lat)...
       +r*K6*(cos(lat))^(R6-1.0)...
       *(R6*sin(lat)*sin(lat)-cos(lat)*cos(lat))*cos(R6*lon); 
    VelSp(2)=-r*K6*R6*(cos(lat))^(R6-1.0)*sin(lat)*sin(R6*lon);
end
end


    