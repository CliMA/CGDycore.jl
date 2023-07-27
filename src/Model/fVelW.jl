function fVelW(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.ProfVelW)
  if str == "linear"
    w=x[1];
  elseif str == "advectionspheredcmip"
    (Lon,Lat,R) = cart2sphere(x[1],x[2],x[3])
    Z=max(R-Phys.RadEarth,0);
    pZ = Phys.p0 * exp(-Z / Param.ScaleHeight)
    RhoZ = pZ / Phys.Rd / Param.T_0
    LonP = Lon - 2 * pi * time / Param.tau
    k = 10 * Phys.RadEarth / Param.tau
    sp =
       1 + exp((Param.p_top - Phys.p0) / Param.b / Param.p_top) - exp((pZ - Phys.p0) / Param.b / Param.p_top) -
         exp((Param.p_top - pZ) / Param.b / Param.p_top)
    omega = Param.omega_0 * sin(LonP) * cos(Lat) * cos(2 * pi * time / Param.tau) * sp
    w = -omega / RhoZ / Phys.Grav  
  elseif str == "advectiontestdeform"
    w = Param.uMax * sinpi(x[3] / Param.H) * cospi(time / Param.EndTime)
  elseif str == "heldsuarezcart"
    w=sin(pi*x[3]/Param.H)
  elseif str == "rotationalcart"
    w = sinpi(x[3] / Param.H) * cospi(time / Param.EndTime)
#   w = (x[1] - Param.xC) / 1000.0
  else
    w=0;
  end
  return w
end

