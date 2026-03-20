function ConstructDG(k,NodalPoints,ElemType::Grids.Tri)
  s = @polyvar x[1:2]

  if k == 0
    println("error: k must be greater than 0")
    return
  end  
  phi = DG.Polynomial_k(k,s)
  DoF = length(phi)
  phiB = Array{Polynomial,1}(undef,DoF)
  Gradphi = Array{Polynomial,2}(undef,DoF,2)
  I = zeros(DoF,DoF)
# Compute functional over nodes
  @inbounds for iDoF = 1 : DoF
    @inbounds for jDoF = 1 : DoF
      I[iDoF,jDoF] = phi[jDoF](NodalPoints[1,iDoF],NodalPoints[2,iDoF])
    end
  end  
  @inbounds for iDoF = 1 : DoF  
    @inbounds for jDoF = 1 : DoF  
      if abs(I[iDoF,jDoF]) < 1.e-12
        I[iDoF,jDoF] = 0
      end
    end
  end  
  r = zeros(DoF)
  @inbounds for iDoF = 1 : DoF  
    r[iDoF] = 1
    c = I \ r
    phiB[iDoF] = 0.0 * x[1] + 0.0 * x[2]
    @inbounds for jDoF = 1 : DoF  
      phiB[iDoF] += c[jDoF] * phi[jDoF]
    end  
    phiB[iDoF] = round.(phiB[iDoF], digits=5)
    r[iDoF] = 0
    Gradphi[iDoF,1] = differentiate(phiB[iDoF],x[1])
    Gradphi[iDoF,2] = differentiate(phiB[iDoF],x[2])
  end  
  return Gradphi
end

