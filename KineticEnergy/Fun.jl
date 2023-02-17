function FunProjectLG!(c,f,X,Fe)
  OPx = size(c,1)
  OPz = size(c,2)
  z = 0.0
  for i = 1 : OPx
    for j = 1 :OPz
      @views z = sum(Fe.IntZLGL2LG[j,:] .* X[i,:,2])  
      c[i,j] = f(X[i,j,1],z)  
    end
  end
end

function FunProject!(c,f,X)
  OPx = size(c,1)
  OPz = size(c,2)
  z = 0.0
  for i = 1 : OPx
    for j = 1 :OPz
      c[i,j] = f(X[i,j,1],X[i,j,2])
    end
  end
end

function uFun(x,z)
  f=sin(pi*x)*sin(pi*z)
end
function DxuFun(x,z)
# f=sin(pi*x)*sin(pi*z)
  Dxf = pi*cos(pi*x)*sin(pi*z)
end
function DzuFun(x,z)
# f=sin(pi*x)*sin(pi*z)
  Dzf = pi*sin(pi*x)*cos(pi*z)
end  

function wFun(x,z)
  f=sin(2*pi*x)*sin(2*pi*z)
end
function DxwFun(x,z)
# f=sin(2*pi*x)*sin(2*pi*z)
  Dxf=2*pi*cos(2*pi*x)*sin(2*pi*z)
end
function DzwFun(x,z)
# f=sin(2*pi*x)*sin(2*pi*z)
  Dzf=2*pi*sin(2*pi*x)*cos(2*pi*z)
end

function RhoFun(x,z)
  f=sin(2*pi*x)*sin(2*pi*z) + 1.0
end
function DxRhoFun(x,z)
# f=sin(2*pi*x)*sin(2*pi*z) + 1.0
  Dxf=2*pi*cos(2*pi*x)*sin(2*pi*z)
end
function DzRhoFun(x,z)
# f=sin(2*pi*x)*sin(2*pi*z) + 1.0
  Dzf=2*pi*sin(2*pi*x)*cos(2*pi*z)
end

function DivFun(x,z)
# Div f
# d/dx (sin(pi*x)*sin(pi*z)*(sin(2*pi*x)*sin(2*pi*z) + 1.0) +
# d/dz (sin(2*pi*x)*sin(2*pi*z)*(sin(2*pi*x)*sin(2*pi*z) + 1.0)

  u = uFun(x,z)
  Dxu = DxuFun(x,z)
  w = wFun(x,z)
  Dzw = DzwFun(x,z)
  Rho = RhoFun(x,z)
  DxRho = DxRhoFun(x,z)
  DzRho = DzRhoFun(x,z)
  return -(Dxu + Dzw) * Rho - u * DxRho - w * DzRho
end

function GradZKin(x,z)
#  Kin = 0.5*(uFun*uFun+wFun*wFun)
  u = uFun(x,z)
  Dzu = -DzuFun(x,z)
  w = wFun(x,z)
  Dzw = -DzwFun(x,z)
  return u * Dzu + w *Dzw
end   

function GradXKin(x,z)
#  Kin = 0.5*(uFun*uFun+wFun*wFun)
  u = uFun(x,z)
  Dxu = -DxuFun(x,z)
  w = wFun(x,z)
  Dxw = -DxwFun(x,z)
  return u * Dxu + w *Dxw
end

function Curly(x,z)
#  Curly = (Dxw,-Dzu)
  Dzu = DzuFun(x,z)
  Dxw = DxwFun(x,z)
  return Dxw-Dzu
end

function AdvuMom(x,z)
#  w*Curly
  w = wFun(x,z)
  Curl = Curly(x,z)
  return w*Curl
end

function AdvwMom(x,z)
# u*Curly
  u = uFun(x,z)
  Curl = Curly(x,z)
  return -u*Curl
end


