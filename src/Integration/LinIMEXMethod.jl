function LinIMEXMethod{FT}() where FT<:AbstractFloat
  nStage=0
  AE = zeros(FT,0,0)
  bE = zeros(FT,0)
  AI = zeros(FT,0,0)
  bI = zeros(FT,0)
  D = zeros(FT,0,0)
  E = zeros(FT,0,0)
  d = zeros(FT,0)
  e = zeros(FT,0)
# AHat=zeros(0,0)
# BHat=zeros(0)
# A=zeros(0,0)
# B=zeros(0)
  name = ""
  JacComp = false
  return LinIMEXMethod{FT}(
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
#   AHat,
#   BHat,
#   A,
#   B,
    JacComp,
  )
end

function LinIMEXMethod{FT}(IMEX::IMEXDirkMethod) where FT<:AbstractFloat
  str = IMEX.name
  nStage = IMEX.nStage
  ahat = IMEX.AE
  bhat = IMEX.bE
  a = IMEX.AI
  b = IMEX.bI
  Gamma = a[2:nStage,2:nStage] - aHat[2:nStage,2:nStage]
  b[2:nStage] = b[2:nStage] - bHat[2:nStage]
  AHat = zeros(nStage,nStage)
  BHat = zeros(nStage)
  A = zeros(nStage,nStage)
  B = zeros(nStage)
  AHat[2:nStage,1:nStage-1] = Gamma \ aHat[2:nStage,1:nStage-1]
  A[2:nStage,2:nStage] = inv(Gamma)
  BHat[1:nStage-1] = bHat[1:nStage-1] - AHat[2:nStage,1:nStage-1]' * b[2:nStage]
  BHat[nStage] = bHat[nStage]
  B[2:nStage] = inv(Gamma)' * b[2:nStage]
  JacComp = true
  return LinIMEXMethod{FT}(
    str,
    nStage,
    AHat,
    BHat,
    A,
    B,
    JacComp,
  )
end

function LinIMEXMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
  @show str
  if str == "IMEXEuler"
    nStage = 2  
    aHat = [0 0
            1.0 0]
    bHat = [1.0 0]        
    a = [0 0
         0.0 1.0]
    b = [0.0 1.0]        
  elseif str == "ARS343"
    nStage = 4
    gamma = 0.4358665215084590
    b1 = -3/2 * gamma^2 + 4 * gamma - 1/4
    b2 =  3/2 * gamma^2 - 5 * gamma + 5/4
    AI = [0 0 0 0;
         0 gamma 0 0;
        0 (1-gamma)/2 gamma 0;
        0 b1 b2 gamma]
    bI = [0; b1; b2; gamma]

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
    bE = [  0; b1; b2; gamma]
  elseif str == "AR2"
    @show str
    nStage = 3
    s2 = sqrt(2)
    a32 = 1/6*(3+2*s2)
    AE = [0     0    0
            2-s2 0 0
            1-a32 a32 0]
    bE = [1/(2*s2); 1/(2*s2); 1-1/s2]        
    
    AI=[0 0 0
       1-1/s2 1-1/s2 0
       1/(2*s2) 1/(2*s2) 1-1/s2]
    bI = [1/(2*s2); 1/(2*s2); 1-1/s2]        
  elseif str == "ARS222"
    nStage = 3
    gammaD = (2 − sqrt(2))/2
    delta = 1 - 1 / (2 * gammaD)
    AE = [0 0 0
        gammaD 0 0
        delta 1 - delta 0]
    bE = [delta; 1 - delta; 0]
    cE = zeros(FT,nStage)
    AI = [0 0 0
        0 gammaD 0
        0 1 - gammaD gammaD]
    bI = [0; 1 - gammaD; gammaD]  
  elseif str == "M1HOMME"
    nStage = 6
    aHat = [  0   0   0   0   0   0
             1/5  0   0   0   0   0
              0  1/5  0   0   0   0
              0   0  1/3  0   0   0
              0   0   0  1/2  0   0
              0   0   0   0   1   0]
    bHat = [  0   0   0   0   1   0]
    a    = [  0       0     0   0   0    0
              0      1/5    0   0   0    0
              0       0    1/5  0   0    0
              0       0     0  1/3  0    0
              0       0     0   0  1/2   0
              5/18   5/18   0   0   0   8/18]
    b    =  [ 5/18   5/18   0   0   0   8/18]  
  end
# Gamma = a[2:nStage,2:nStage] - aHat[2:nStage,2:nStage]
# b[2:nStage] = b[2:nStage] - bHat[2:nStage]
# AHat = zeros(nStage,nStage)
# BHat = zeros(nStage)
# A = zeros(nStage,nStage)
# B = zeros(nStage)
# AHat[2:nStage,1:nStage-1] = Gamma \ aHat[2:nStage,1:nStage-1]
# A[2:nStage,2:nStage] = inv(Gamma)
# BHat[1:nStage-1] = bHat[1:nStage-1] - AHat[2:nStage,1:nStage-1]' * b[2:nStage]
# BHat[nStage] = bHat[nStage]
# B[2:nStage] = inv(Gamma)' * b[2:nStage]

  AI = AI - AE
  bI = bI - bE
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
  return LinIMEXMethod{FT}(
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
#   AHat,
#   BHat,
#   A,
#   B,
    JacComp,
  )
end
