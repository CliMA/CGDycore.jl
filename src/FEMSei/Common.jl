function CGLine(k,x)
  phi = Array{Polynomial,1}(undef,k+1)
# phi[1] = 1.0 + 0.0 * x
  for i = 0 : k
    phi[i+1] = 0.5^k * (1.0-x)^(k-i)*(1.0+x)^i
#   phi[i+1] = phi[i] * x
  end  
  return phi
end

function Polynomial_1D(k,x,l)
  phi = Array{Polynomial,1}(undef,k+1)
  phi[1] = 0.0 * x[1][1] + 0.0 * x[1][2] + 1.0 
  for i = 1 : k
    phi[i+1] = x[1][l] * phi[i]
  end  
  return phi
end

function HomegenuousPolynomial(k,x)
  phi = Array{Polynomial,1}(undef,k+1)
  for i = 0 : k
    phi[i+1] = x[1][1]^(k-i) * x[1][2]^i + 0.0 
  end  
  return phi
end  

function Polynomial_k(k,x)
  DoF::Int = (k + 1) *(k + 2) / 2
  phi = Array{Polynomial,1}(undef,DoF)
  iDoF = 1
  for i = 0 : k
    for j = 0 : i  
      phi[iDoF] = x[1][1]^(i-j) * x[1][2]^j + 0.0 
      iDoF += 1
    end
  end
  return phi
end

