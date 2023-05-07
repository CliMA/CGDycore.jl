# Note we compute -Grad
function fScalar(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  # global Omega uM lat0G lat1G eN
  str = lowercase(Model.ProfTest)
  if str == "testgrad"  
    Lx = 10000.0  
    Lz = 10000.0
    c = (x[1]/Lx)^3+(x[3]/Lz)^3
  elseif str == "testgrad1"  
    Lx = 10000.0
    Lz = 10000.0
    c = cos(pi * x[1] / Lx) * cos(pi * x[3] / Lz)
  end
  return c
end
function fGrad12(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  # global Omega uM lat0G lat1G eN
  str = lowercase(Model.ProfTest)
  if str == "testgrad"  
    Lx = 10000.0  
    Lz = 10000.0
    dcdx1 = -3.0 * (x[1]/Lx)^2 / Lx
    dcdx2 = -0.0
  elseif str == "testgrad1"  
    Lx = 10000.0
    Lz = 10000.0
#   c = cos(pi * x[1] / Lx) * cos(pi * x[3] / Lz)
    dcdx1 = sin(pi * x[1] / Lx) *pi /Lx  * cos(pi * x[3] / Lz)
    dcdx2 = -0.0
  elseif str == "testgrad2"  
    Lx = 10000.0
    Lz = 10000.0
    c = cos(pi * x[1] / Lx) * cos(pi * x[3] / Lz)
    dcdx1 = sin(pi * x[1] / Lx) *pi /Lx  * cos(pi * x[3] / Lz)
    dcdx1 = c * dcdx1
    dcdx2 = 0.0
  end
  return dcdx1,dcdx2
end  

function fGrad3(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  # global Omega uM lat0G lat1G eN
  str = lowercase(Model.ProfTest)
  if str == "testgrad"
    Lx = 10000.0  
    Lz = 10000.0
    dcdx3 = -3.0 * (x[3]/Lz)^2 / Lz
  elseif str == "testgrad1"  
    Lx = 10000.0
    Lz = 10000.0
#   c = cos(pi * x[1] / Lx) * cos(pi * x[3] / Lz)
    dcdx3 = cos(pi * x[1] / Lx) * sin(pi * x[3] / Lz) * pi / Lz
  elseif str == "testgrad2"  
    Lx = 10000.0
    Lz = 10000.0
    c = cos(pi * x[1] / Lx) * cos(pi * x[3] / Lz)
    dcdx3 = cos(pi * x[1] / Lx) * sin(pi * x[3] / Lz) * pi / Lz
    dcdx3 = c * dcdx3
  end
  return dcdx3
end

