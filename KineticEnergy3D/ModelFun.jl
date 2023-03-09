function fProjectC!(c,f,X,Fe,Phys,Param)
  OPx = size(c,1)
  OPy = size(c,2)
  for i = 1 : OPx
    for j = 1 : OPy
      x  = 0.5 * (X[i,j,1,1] + X[i,j,2,1])
      y  = 0.5 * (X[i,j,1,2] + X[i,j,2,2])
      z  = 0.5 * (X[i,j,1,3] + X[i,j,2,3])
      c[i,j] = f(x,y,z,Phys,Param)
    end
  end
end

function fRho(x,y,z,Phys,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    Grav =Phys.Grav
    p0 = Phys.p0
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    DeltaTh = Param.DeltaTh
    xC0 = Param.xC0
    zC0 = Param.zC0
    rC0 = Param.rC0
    pLoc =p0*(1-Grav*z*kappa/(Rd*Th0))^(1/kappa)
    rr = sqrt((x-xC0)^2+(z-zC0)^2)
    ThLoc = Th0
    if rr <rC0
      ThLoc =ThLoc+DeltaTh*cos(0.5*pi*rr/rC0)^2
    end
    Rho = pLoc / ((pLoc / p0)^kappa * Rd * ThLoc)
  elseif Example == "AdvectionCart"
    Rho = 1.0
  elseif Example == "HillAgnesiXCart" || Example == "HillAgnesiYCart" 
    NBr = Param.NBr
    Grav = Phys.Grav
    p0 = Phys.p0
    Cpd = Phys.Cpd
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    S = NBr * NBr / Grav
    ThLoc = Th0 * exp(z * S)
    pLoc = p0 *(1.0 - Grav/(Cpd * Th0 * S) * (1.0 - exp(-S * z)))^(Cpd / Rd)
    Rho = pLoc / ((pLoc / p0).^kappa * Rd*ThLoc)
  end    
  return Rho
end

function fTheta(x,y,z,Phys,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    Grav =Phys.Grav
    p0 = Phys.p0
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    DeltaTh = Param.DeltaTh
    xC0 = Param.xC0
    zC0 = Param.zC0
    rC0 = Param.rC0
    pLoc =p0*(1-Grav*z*kappa/(Rd*Th0))^(1/kappa)
    rr = sqrt((x-xC0)^2+(z-zC0)^2)
    ThLoc = Th0
    if rr <rC0
      ThLoc =ThLoc+DeltaTh*cos(0.5*pi*rr/rC0)^2
    end
    Theta = ThLoc
  elseif Example == "HillAgnesiXCart" || Example == "HillAgnesiYCart" 
    NBr = Param.NBr
    Grav = Phys.Grav
    p0 = Phys.p0
    Cpd = Phys.Cpd
    Rd = Phys.Rd
    kappa = Phys.kappa
    Th0 = Param.Th0
    S = NBr * NBr / Grav
    ThLoc = Th0 * exp(z * S)
    Theta = ThLoc
  elseif Example == "AdvectionCart"
    r = sqrt((x - 0.25)^2 + (y - 0.25)^2)
    if r <= 0.1 
      Theta = 1.0
    else
      Theta = 0.0
    end  
  end    
  return Theta
end

function fuVel(x,y,z,Phys,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    u = Param.uMax  
  elseif Example == "HillAgnesiXCart" || Example == "HillAgnesiYCart" 
    u = Param.uMax  
  elseif Example == "AdvectionCart"
    u = -y 
  end
  return u
end

function fvVel(x,y,z,Phys,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    v = Param.vMax  
  elseif Example == "HillAgnesiXCart" || Example == "HillAgnesiYCart" 
    v = Param.vMax  
  elseif Example == "AdvectionCart"
    v = x
  end
  return v
end

function fwVel(x,y,z,Phys,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    w = Param.wMax  
  elseif Example == "HillAgnesiXCart" || Example == "HillAgnesiYCart" 
    w = Param.wMax  
  end
  return w
end


