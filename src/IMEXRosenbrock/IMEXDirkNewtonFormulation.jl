function IMEXDirkNewtonFormulation(IMEX)

# y = y_n + h*AE*fE + h*AI*fI
# or
# Z = h*AE*fE + h*AI*fI
# Z_K = (Z_2,Z_3,...,Z_k)
# fE_K=(fE_1,...,fE_{k-1})
# fI_K=(fE_2,...,fE_k})

# Stage 1 to i
# Z_K = h*AE_K*fE_K + h*AI_K*fI_K
# AI_K^{-1}*Z_K-h*AI_K^{-1}*AE_K*fE_K = fI_K

# Z_{k+1} = AE[k+1,1:k]*fE_K + AI[k+1,2:k]*fI_K + fI(yn + Z_{k+1})
# Z_{k+1} = AAE[k+1,1:k]*fE_K + AII[k+1,2:k]*Z_K + fI(yn + Z_{k+1})
# AAE[k+1,1:k] = AE[k+1,1:k] - AI[k+1,2:k] * AI_K^{-1}*AE_K
# AII[k+1,1:k] = AI_K^{-1}{-1}


  FT = eltype(IMEX.AE)

  nStage = IMEX.nStage - 1

  AE = [IMEX.AE
        IMEX.bE']
  AI = [IMEX.AI
        IMEX.bI']      

  AAE = zeros(size(AE))      
  AAI = zeros(size(AI))      

  AAE[2,:] = AE[2,:]
  AAI[2,:] = AI[2,:]


  @show AI

  for k = 2 : IMEX.nStage
    AAE[k+1:k+1,1:k] = AE[k+1:k+1,1:k] - AI[k+1:k+1,2:k] * inv(AI[2:k,2:k]) * AE[2:k,1:k]  
    AAI[k+1:k+1,2:k] = AI[k+1:k+1,2:k] * inv(AI[2:k,2:k])   
  end    
  return AAE,AAI
end
