function IMEXDirkToRosenbrock(IMEX)

  FT = eltype(IMEX.AE)

  nStage = IMEX.nStage - 1 

  @views alphaR = IMEX.AE[2:end,1:end-1]
  @views d = IMEX.AI[2:end,2:end] - IMEX.AE[2:end,2:end]
  @views b = IMEX.bE[1:end-1]

  gamma = inv(alphaR) * d * alphaR
  alpha = IMEX.AE[1:end-1,1:end-1]
  gammaD = gamma[2,2]

  a = alpha / gamma
  c = -inv(gamma)
  @show size(gamma)
  @show size(b)
  m = gamma'\b

  return INT.RosenbrockMethod{FT}(
    "ROS"*IMEX.name,
    nStage,
    # Standard formulation
    alpha,
    gamma,
    b,  
    # jacobian free formulation
    a,
    c,
    gammaD,
    m,
  )  
end    
