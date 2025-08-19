using LinearAlgebra
function LocalInterpolation(ILoc,F)
  @. ILoc = 0
  grad = zeros(2)
  ILoc[1,1] = 1
  n = [1;0];
  ksi1 = 1.0
  ksi2 = -1.0
  _,J,dXdxIT = Jacobi(F,ksi1,ksi2)
  e = norm(F.P2-F.P3)
  for i = 1 : 4
    gradphi(grad,ksi1,ksi2,i)
    ILoc[2,i] = 2.0 * J / e * (dXdxIT * n)' * (dXdxIT * grad)
  end

  n =[0;1]
  ksi1 = -1.0
  ksi2 = 1.0
  _,J,dXdxIT = Jacobi(F,ksi1,ksi2)
  e = norm(F.P4-F.P3)
  for i = 1 : 4
    gradphi(grad,ksi1,ksi2,i)
    ILoc[4,i] = 2.0 * J / e * (dXdxIT * n)' * (dXdxIT * grad)
  end
end

function phi(ksi1,ksi2,i)
  if i == 1
    f = 1/4*(1-ksi1)*(1-ksi2)
  elseif i == 2  
    f = 1/4*(1+ksi1)*(1-ksi2)
  elseif i == 3
    f = 1/4*(1+ksi1)*(1+ksi2)
  elseif i == 4
    f = 1/4*(1-ksi1)*(1+ksi2)
  end
  return f
end

function gradphi(f,ksi1,ksi2,i)
  if i == 1
    f[1,1] = -1/4*(1-ksi2)
    f[2,1] = -1/4*(1-ksi1)
  elseif i == 2  
    f[1,1] = 1/4*(1-ksi2)
    f[2,1] = -1/4*(1+ksi1)
  elseif i == 3
    f[1,1] = 1/4*(1+ksi2)
    f[2,1] = 1/4*(1+ksi1)
  elseif i == 4
    f[1,1] = -1/4*(1+ksi2)
    f[2,1] = 1/4*(1-ksi1)
  end
end

function Jacobi(F,ksi1,ksi2)
  dXdx = zeros(2,2)
  dXdx[1,1] = 0.25*(-F.P1[1]*(1-ksi2)+F.P2[1]*(1-ksi2)+F.P3[1]*(1+ksi2)-F.P4[1]*(1+ksi2))
  dXdx[2,1] = 0.25*(-F.P1[2]*(1-ksi2)+F.P2[2]*(1-ksi2)+F.P3[2]*(1+ksi2)-F.P4[2]*(1+ksi2))

  dXdx[1,2] = 0.25*(-F.P1[1]*(1-ksi1)-F.P2[1]*(1+ksi1)+F.P3[1]*(1+ksi1)+F.P4[1]*(1-ksi1))
  dXdx[2,2] = 0.25*(-F.P1[2]*(1-ksi1)-F.P2[2]*(1+ksi1)+F.P3[2]*(1+ksi1)+F.P4[2]*(1-ksi1))
  J = det(dXdx)
  dXdxIT = inv(dXdx')
  return (dXdx,J,dXdxIT)
end

pert = 0.1
pert = 0.0
PC = [0+pert,0+pert]
P1 = [1+pert,1+pert]
P2 = [-1,1]
P3 = [-1,-1]
P4 = [1,-1]

PM1 = [0,1]
PM2 = [-1,0]
PM3 = [0,-1]
PM4 = [1,0]
mutable struct FaceT
  P1::Array{Float64, 1}
  P2::Array{Float64, 1}
  P3::Array{Float64, 1}
  P4::Array{Float64, 1}
end  
F = Array{FaceT,1}(undef,4)
F[1]=FaceT(P1,PM1,PC,PM4)
F[2]=FaceT(P2,PM2,PC,PM1)
F[3]=FaceT(P3,PM3,PC,PM2)
F[4]=FaceT(P4,PM4,PC,PM3)
ILoc = zeros(4,4,4)
GlobLoc = zeros(Int,4,4)
ksi1 = 1
ksi2 = -1
I = zeros(9,9)
for i = 1 : 4
  GlobLoc[:,i] = [i,4+i,9,4+i-1]
  if i == 1
    GlobLoc[4,i] = 8
  end
  @views LocalInterpolation(ILoc[:,:,i],F[i])
  @show ILoc[:,:,i]
  @. I[GlobLoc[:,i],GlobLoc[:,i]] += ILoc[:,:,i]
end
NumF = 4
Grad = zeros(NumF,NumF)
I[2*NumF+1,2*NumF+1] = -4
I[2*NumF+1,1] = 1
I[2*NumF+1,2] = 1
I[2*NumF+1,3] = 1
I[2*NumF+1,4] = 1
for i = 1 : NumF
  b = zeros(2*NumF+1)
  b[i] = 1  
  b[2*NumF+1] = 0.0
# IR = I[1:2*NumF,:]
# c = IR' * ((IR*IR') \ b[1:2*NumF])
  c = I \ b
  for j = 1 : NumF
     Grad[j,i] = sum(ILoc[2,:,j] .* c[GlobLoc[:,j]])
  end
end
nothing


