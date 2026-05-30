@inline function ComputeVolume(J,wH,wV)
  DoF = size(J,1)
  M = size(J,2)
  Vol = 0.0
  @inbounds for iDoF = 1 : DoF
    Vol3 = 0.0  
    @inbounds for k = 1 : M
      Vol3 += wV[k] / J[iDoF,k]
    end
    Vol += wH[iDoF] * Vol3
  end
  return Vol
end  

@inline function ComputeLowerSurface(dXdxI,wH)
  DoF = size(dXdxI,2)
  Area = 0.0
  @inbounds for iDoF = 1 : DoF
    @views s = norm(dXdxI[:,iDoF])  
    Area += wH[iDoF] * s 
  end
  return Area 
end

@inline function ComputeSideXSurface(dXdxI,wH,wV,I)
  M = size(wV,1)
  N = size(wH,1)
  Area = 0.0
  @inbounds for j = 1 : N
    iDoF = I + (j - 1) * N
    Area1 = 0.0
    @show iDoF
    @inbounds for k = 1 : M
      @views s = norm(dXdxI[:,k,iDoF])
      Area1 += wV[k] * s
    end  
    Area += wH[j] * Area1
  end
  return Area
end

@inline function ComputeSideYSurface(dXdxI,wH,wV,J)
  M = size(wV,1)
  N = size(wH,1)
  Area = 0.0
  @inbounds for i = 1 : N
    iDoF = i + (J - 1) * N
    Area1 = 0.0
    @inbounds for k = 1 : M
      @views s = norm(dXdxI[:,k,iDoF])
      Area1 += wV[k] * s
    end
    Area += wH[i] * Area1
  end
  return Area
end

