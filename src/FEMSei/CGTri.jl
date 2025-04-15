#constructing the CG elements for the Triangular grid
function ConstructCG(k,ElemType::Grids.Tri)
  s = @polyvar x[1:2]

  if k == 0
    stop
  end  
  P_k = DG.Polynomial_k(k,s)
  lP_k = length(P_k)
  DoF = lP_k 
  DoFE = k - 1
  DoFN = 1
  DoFF = DoF - 3 * DoFE - 3 * DoFN 
  phi = Array{Polynomial,1}(undef,DoF)
  phiB = Array{Polynomial,1}(undef,DoF)
  iDoF = 1 
  @inbounds for i = 1 : lP_k
    phi[i] = P_k[i]  
  end  
  @polyvar t
  I = zeros(DoF,DoF)
  rDoF = 1
# Compute functional over nodes
  @inbounds for iDoF = 1 : DoF
    I[rDoF,iDoF] = phi[iDoF](-1.0,-1.0)
    I[rDoF+1,iDoF] = phi[iDoF](1.0,-1.0)
    I[rDoF+2,iDoF] = phi[iDoF](-1.0,1.0)
  end
  rDoF += 3
# Compute functional over Edges
  # Edge 1 (-1,-1) -> (1,-1)
  @inbounds for iDoF = 1 : DoF
    @inbounds for i = 0 : k - 2
      ksi1 = -1.0 + 2 * (i + 1) / k
      ksi2 = -1.0
      if iDoF == 1
        @show "E1",ksi1,ksi2
      end
      I[rDoF+i,iDoF] = phi[iDoF](ksi1,ksi2)
    end
  end
  rDoF += k - 1
  # Edge 2 (1,-1) -> (-1,1)
  @inbounds for iDoF = 1 : DoF
    @inbounds for i = 0 : k - 2
      ksi1 = 1.0 - 2 * (i + 1) / k
      ksi2 = -1.0 + 2 * (i + 1) / k
      I[rDoF+i,iDoF] = phi[iDoF](ksi1,ksi2)
    end
  end
  rDoF += k - 1
  # Edge 3 (-1,-1) -> (1,1)
  @inbounds for iDoF = 1 : DoF
    @inbounds for i = 0 : k - 2
      ksi1 = -1.0 
      ksi2 = -1.0 + 2 * (i + 1) / k
      I[rDoF+i,iDoF] = phi[iDoF](ksi1,ksi2)
    end
  end
  rDoF += k - 1
  @show DoF,DoFF
  @show rDoF
# Interior  
  @inbounds for iDoF = 1 : DoF
    incr = 0
    @inbounds for j = 0 : k - 2
      @inbounds for i = 1  : k - 2 - j
        # Lage der Punkte
        ksi1 = -1.0 + 2 * i / k 
        ksi2 = -1.0 + 2 * (j + 1) / k
        if iDoF == 1
          @show ksi1,ksi2
        end  
        I[rDoF+incr,iDoF] = phi[iDoF](ksi1,ksi2)
        incr += 1
      end  
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
  end  
  return DoF, DoFE, DoFF, phiB
end

