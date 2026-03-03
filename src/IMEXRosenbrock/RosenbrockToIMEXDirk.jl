function RosenbrockToIMEXDirk(ROS)

  FT = eltype(ROS.alpha)

  # Y1 = y
  # k1 = h*f(Y1) + h*J(gamma11*k1)
  # Y2 = y + alpha21*k1
  # k2 = h*f(Y2) + h*J(gamma21*k1+gamma22*k2)
  # Y3 = y + alpha31*k1 + alpha32*k2
  # k3 = h*f(Y3) + h*J(gamma31*k1+gamma32*k2+gamma33*k3)
  # Y4 = y + alpha41*k1 + alpha42*k2 + alpha43*k3
  # k4 = h*f(Y4) + h*J(gamma41*k1+gamma42*k2+gamma43*k3+gamma44*k4)
  # yNew = b1*k1 + b2*k2 + b3*k3 + b4 * k4

  nStage = ROS.nStage + 1
  alphaR = [ROS.alpha[2:end,1:end]
            ROS.b']
  AE = [zeros(FT,1,nStage)
        alphaR zeros(FT,nStage-1)]
  bE = AE[end,:]

  AI = [zeros(FT,1,nStage)
        zeros(FT,nStage-1) alphaR * ROS.gamma * inv(alphaR)]
  AI[:,2:end] .+= AE[:,2:end]      
  bI = AI[end,:]
  cE = zeros(FT,nStage)
  cI = zeros(FT,nStage)
  for i = 1 : nStage
    cE[i] = sum(AE[i,:])  
    cI[i] = cE[i]
  end  
  return IMEXDirkMethod{FT}(
    "IMEX"*ROS.name,
    nStage,
    AE,
    bE,
    cE,
    AI,
    bI,
    cI,
  )
end
