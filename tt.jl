using DynamicPolynomials
using FastGaussQuadrature
function Lagrange(x,xw,i)
  L = 1 + 0*x
  for j = 1 : length(xw)
    if j != i
      L = L * (x - xw[j]) / (xw[i] - xw[j])
    end
  end
  return L
end

@polyvar x[1:2] 
k = 3
L1 = Array{Polynomial,1}(undef,k)
L2 = Array{Polynomial,1}(undef,k)
xw,_= gausslobatto(3)
for i = 1 : k
  L1[i] = Lagrange(x[1],xw,i)  
  L2[i] = Lagrange(x[2],xw,i)  
end  

DoF = k * k
phi = Array{Polynomial,2}(undef,DoF,1)
for j = 1 : k
  for i = 1 : k
    phi[iDoF,1] = L1[i] * L2[j]  
    iDoF += 1
  end
end  
    
