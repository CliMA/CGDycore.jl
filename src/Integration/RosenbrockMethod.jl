function RosenbrockMethod{FT}() where FT<:AbstractFloat
  nStage = 0
  alpha = zeros(FT,0,0)
  gamma = zeros(FT,0,0)
  b = zeros(FT,0)
  gammaD = FT(0)
  a = zeros(FT,0,0)
  c = zeros(FT,0,0)
  m = zeros(FT,0)
  name = ""
  JacComp = true
  return RosenbrockMethod{FT}(
    name,
    nStage,
    alpha,
    gamma,
    b,
    a,
    c,
    gammaD,
    m,
    JacComp,
  )
end

function RosenbrockMethod{FT}(RK::RungeKuttaExMethod,gammaD,gammaV) where FT<:AbstractFloat
  nStage = RK.nStage
  alpha = RK.A
  b = RK.b
  iV = 1
  gamma = zeros(nStage,nStage)
  for iStage = 1 : nStage
    gamma[iStage,iStage] = gammaD
    for jStage = 1 : iStage - 1
       gamma[iStage,jStage] = gammaV[iV]
       iV += 1
    end
  end
  a = alpha / gamma
  c = -inv(gamma)
  m = gamma'\b
  JacComp = true
  return RosenbrockMethod{FT}(
    "ROS"*RK.name,
    nStage,
    alpha,
    gamma,
    b,
    a,
    c,
    gammaD,
    m,
    JacComp,
  )
end


function RosenbrockMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
  JacComp = true
  if str == "SSP-Knoth"
    nStage = 3
    alpha = zeros(FT,nStage,nStage)
    alpha[2,1] = 1
    alpha[3,1] = 1/4
    alpha[3,2] = 1/4
    b = zeros(FT,nStage)
    b[1] = 1/6
    b[2] = 1/6
    b[3] = 2/3
    Gamma = zeros(FT,nStage,nStage)
    Gamma[1,1] = 1
    Gamma[2,2] = 1
    Gamma[3,1] = -3/4
    Gamma[3,2] = -3/4
    Gamma[3,3] = 1
    gamma = FT(1)        
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']
  elseif str == "ROS2W"
    @show "ROS2W"
    nStage = 2
    gamma = 1 + 1 / sqrt(2)
    alpha = zeros(FT,nStage,nStage)
    b = zeros(FT,nStage)
    b[2] = 0.5
    b[1] = 1 - b[2]
    alpha[2,1] = 1 / (2 * b[2])
    Gamma = zeros(FT,nStage,nStage)
    Gamma[1,1] = gamma
    Gamma[2,1] = -b[2] / gamma
    Gamma[2,2] = gamma
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']
  elseif str == "RosEul"
    nStage = 1
    alpha = zeros(FT,nStage,nStage)
    b = zeros(FT,nStage)
    b[1] = 1
    Gamma = zeros(FT,nStage,nStage)
    Gamma[1,1] = 1
    gamma = FT(1)        
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']
  elseif str == "ROSRK3"
    nStage = 3
    alpha = zeros(FT,nStage,nStage)
    alpha[2,1] = 1/3
    alpha[3,2] = 1/2
    b = zeros(FT,nStage)
    b[3] = 1
    gamma = 1
    Gamma = zeros(FT,nStage,nStage)
    Gamma[1,1] = gamma
    Gamma[2,1] = (1-12*gamma^2) /(-9+36*gamma)
    Gamma[2,2] = gamma
    Gamma[3,1] = -1/4+2*gamma
    Gamma[3,2] = 1/4-3*gamma
    Gamma[3,3] = gamma
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']
  elseif str == "ARS343"
    nStage = 3
    gamma = 0.4358665215084590
    a42 = 0.5529291480359398
    a43 = a42
    b1 = -3 / 2 * gamma^2 + 4 * gamma - 1 / 4
    b2 = 3 / 2 * gamma^2 - 5 * gamma + 5 / 4
    a31 =
        (1 - 9 / 2 * gamma + 3 / 2 * gamma^2) * a42 + 
        (11 / 4 - 21 / 2 * gamma + 15 / 4 * gamma^2) * a43 - 7 / 2 + 13 * gamma - 9 / 2 * gamma^2
    a32 =
        (-1 + 9 / 2 * gamma - 3 / 2 * gamma^2) * a42 + 
        (-11 / 4 + 21 / 2 * gamma - 15 / 4 * gamma^2) * a43 + 4 - 25 / 2 * gamma +
        9 / 2 * gamma^2
    a41 = 1 - a42 - a43
    AHat = [0 0 0 0
            gamma 0 0 0
            a31 a32 0 0
            a41 a42 a43 0]
    bHat = [0, b1, b2, gamma]
    A = [0 0 0 0
         0 gamma 0 0
         0 (1 - gamma)/2 gamma 0
         0 b1 b2 gamma]
    b = [0, b1, b2, gamma]
    alpha = zeros(FT,nStage,nStage)
    @views @. alpha = AHat[2:end,1:end-1]
    Gamma = zeros(FT,nStage,nStage)
    @views @. Gamma = A[2:end,2:end] - AHat[2:end,2:end]

    Gamma = alpha \ Gamma * alpha
    alpha = AHat[1:end,1:end-1]
    a = alpha / Gamma
    c = -inv(Gamma)
    m = a[end,:]
  elseif str == "M1HOMME"
    nStage = 5
    AHat = [  0   0   0   0   0   0
             1/5  0   0   0   0   0
              0  1/5  0   0   0   0
              0   0  1/3  0   0   0
              0   0   0  1/2  0   0
              0   0   0   0   1   0]
    bHat = [  0   0   0   0   1   0]          
    A    = [  0       0     0   0   0    0
              0      1/5    0   0   0    0
              0       0    1/5  0   0    0
              0       0     0  1/3  0    0
              0       0     0   0  1/2   0
              5/18   5/18   0   0   0   8/18]
    b    =  [ 5/18   5/18   0   0   0   8/18]
    alpha = KernelAbstractions.zeros(backend,FT,nStage,nStage)
    @views @. alpha = AHat[2:end,1:end-1]
    Gamma = KernelAbstractions.zeros(backend,FT,nStage,nStage)
    @views @. Gamma = A[2:end,2:end] - AHat[2:end,2:end]
    b = KernelAbstractions.zeros(backend,FT,nStage)
    b[nStage-1] = [0,0,0,0,1,0]          

    Gamma = alpha \ Gamma * alpha       
    alpha = AHat[1:end,1:end-1]

    a = alpha / Gamma
    c = -inv(Gamma)
    m = a[end,:]
	
  elseif str == "SSP-MaExperimental"
    nStage = 4
    alpha = zeros(FT,nStage,nStage)
    alpha[3,2] = 1
    alpha[4,2] = 1/4
    alpha[4,3] = 1/4
    b = zeros(FT,nStage)
    b[2] = 1/6
    b[3] = 1/6
    b[4] = 2/3
    Gamma = zeros(FT,nStage,nStage)
    Gamma[1,1] = 1
    Gamma[2,2] = 1
    Gamma[2,1] = -sqrt(2)
    Gamma[3,2] = 2*sqrt(2) - 3
    Gamma[3,3] = 1
    Gamma[4,1] = (3 - sqrt(2))/4
    Gamma[4,2] = -3/4
    Gamma[4,3] = -3/4
    Gamma[4,4] = 1
    gamma = FT(1)        
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']
  end
  return RosenbrockMethod{FT}(
    str,
    nStage,
    alpha,
    Gamma,
    b,
    a,
    c,
    gamma,
    m,
    JacComp,
  )
end
