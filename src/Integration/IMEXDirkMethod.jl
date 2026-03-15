function IMEXDirkMethod{FT}() where FT<:AbstractFloat
  name = ""
  nStage = 0
  AE = zeros(FT,0,0)
  bE = zeros(FT,0)
  AI = zeros(FT,0,0)
  bI = zeros(FT,0)
  D = zeros(FT,0,0)
  E = zeros(FT,0,0)
  d = zeros(FT,0)
  e = zeros(FT,0)
  JacComp = false
  return IMEXDirkMethod{FT}(
    name,
    nStage,
    AE,
    bE,
    AI,
    bI,
    D,
    E,
    d,
    e,
    JacComp,
  )
end

function IMEXDirkMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
   
  if str == "ExImpEul"
    nStage = 2
    AI = [0 0 ;
         0 1.0]
    bI = [0, 1.0]
    AE = [0 0 ;
         1.0 0]
    bE = [0, 1.0]
  elseif str == "ARS343"
    nStage = 4
    gamma = 0.4358665215084590
    b1 = -3/2 * gamma^2 + 4 * gamma - 1/4
    b2 =  3/2 * gamma^2 - 5 * gamma + 5/4
    AI = [0 0 0 0;
         0 gamma 0 0;
        0 (1-gamma)/2 gamma 0;
        0 b1 b2 gamma]
    bI = [0, b1, b2, gamma]

    dA42 = 0.5529291480359398;
    dA43 = 0.5529291480359398;
    dA31 =(
        (1.0 - 4.5 * gamma + 1.5 * gamma * gamma) * dA42
          + (2.75 - 10.5 * gamma + 3.75 * gamma * gamma) * dA43
          - 3.5 + 13 * gamma - 4.5 * gamma * gamma)
    dA32 = (
        (-1.0 + 4.5 * gamma - 1.5 * gamma * gamma) * dA42
          + (-2.75 + 10.5 * gamma - 3.75 * gamma * gamma) * dA43
          + 4.0 - 12.5 * gamma + 4.5 * gamma * gamma)
    dA41 = 1.0 - dA42 - dA43;
    # explicit tableau
    AE = [
      0 0 0 0;
      gamma 0 0 0;
      dA31 dA32 0 0;
      dA41 dA42 dA43 0]
    bE = [  0, b1, b2, gamma]
    D=zeros(nStage,nStage)
    E=zeros(nStage,nStage)
    d=zeros(nStage)
    e=zeros(nStage)
    D[2,1] = AE[2,1] 

    D[3,1] = AE[3,1] - D[2,1] * AI[3,2] / AI[2,2]
    D[3,2] = AE[3,2] 
    E[3,2] = AI[3,2] / AI[2,2]

    D[4,1] = AE[4,1] - D[2,1] * AI[4,2] / AI[2,2] - D[3,1] * AI[4,3] / AI[3,3]
    D[4,2] = AE[4,2]                              - D[3,2] * AI[4,3] / AI[3,3]
    D[4,3] = AE[4,3] 
    E[4,2] = AI[4,2] / AI[2,2] - E[3,2] * AI[4,3] / AI[3,3]
    E[4,3] = AI[4,3] / AI[3,3]

    d[1] = bE[1] - D[2,1] * bI[2] / AI[2,2] - D[3,1] * bI[3] / AI[3,3] - D[4,1] * bI[4] / AI[4,4]
    d[2] = bE[2]                            - D[3,2] * bI[3] / AI[3,3] - D[4,2] * bI[4] / AI[4,4]
    d[3] = bE[3]                                                       - D[4,3] * bI[4] / AI[4,4]
    d[4] = bE[4] 
    e[2] = bI[2] / AI[2,2] - bI[3] * E[3,2] / AI[3,3]  - bI[4] * E[4,2] / AI[4,4]
    e[3] = bI[3] / AI[3,3]                            - bI[4] * E[4,3] / AI[4,4]
    e[4] = bI[4] / AI[4,4]
  elseif str == "AR2"
    nStage = 3
    s2 = sqrt(2)
    a32 = 1/6*(3+2*s2)
    AE = [0     0    0
          2-s2 0 0
          1-a32 a32 0]
    bE = [1/(2*s2) 1/(2*s2) 1-1/s2]        
    
    AI=[0 0 0
       1-1/s2 1-1/s2 0
       1/(2*s2) 1/(2*s2) 1-1/s2]
    bI = [1/(2*s2) 1/(2*s2) 1-1/s2]        
  end
  AAE = [AE
        bE']
  AAI = [AI
        bI']

  AAAE = zeros(size(AAE))
  AAAI = zeros(size(AAI))

  AAAE[2,:] = AAE[2,:]

  for k = 2 : nStage
    AAAE[k+1:k+1,1:k] = AAE[k+1:k+1,1:k] - AAI[k+1:k+1,2:k] * inv(AAI[2:k,2:k]) * AAE[2:k,1:k]
    AAAI[k+1:k+1,2:k] = AAI[k+1:k+1,2:k] * inv(AAI[2:k,2:k])
  end
  D = AAAE[1:end-1,:]
  d = AAAE[end,:]
  E = AAAI[1:end-1,:]
  e = AAAI[end,:]
  JacComp = true
  return IMEXDirkMethod{FT}(
    str,
    nStage,
    AE,
    bE,
    AI,
    bI,
    D,
    E,
    d,
    e,
    JacComp,
  )
end
