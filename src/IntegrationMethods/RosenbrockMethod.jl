mutable struct RosenbrockStruct
  nStage::Int
  alpha::Array{Float64, 2}
  Gamma::Array{Float64, 2}
  b::Array{Float64, 1}
  a::Array{Float64, 2}
  c::Array{Float64, 2}
  gamma::Float64
  m::Array{Float64, 1}
end
function RosenbrockMethod()
  nStage = 0
  alpha = zeros(0,0)
  Gamma = zeros(0,0)
  b = zeros(0)
  a = zeros(0,0)
  c =zeros(0,0)
  gamma = 0
  m = zeros(0)
  d = zeros(0)
  SSP = SSPRungeKuttaMethod()
  return RosenbrockStruct(
    nStage,
    alpha,
    Gamma,
    b,
    a,
    c,
    gamma,
    m,
  )
end


function RosenbrockMethod(Method)
  str = Method
  if str == "SSP-Knoth"
    nStage = 3
    alpha = [ 0  0  0
              1  0  0
            1/4 1/4 0]
    b = [1/6,1/6,2/3]
    Gamma = [ 1    0  0
              0    1  0
            -3/4 -3/4 1]
    gamma = 1        
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']

  elseif str == "ROSRK3"
    nStage = 3
    alpha=[ 0  0  0
           1/3 0  0
            0 1/2 0]
    gamma = 1
    Gamma=[                        gamma           0     0
           (1-12*gamma^2) /(-9+36*gamma)       gamma     0
                            -1/4+2*gamma 1/4-3*gamma gamma]
    b=[0,0,1]
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']
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

    alpha = AHat[2:end,1:end-1]
    Gamma = A[2:end,2:end] - AHat[2:end,2:end]

    b = [0,0,0,0,1,0]          

    Gamma = alpha \ Gamma * alpha       
    alpha = AHat[1:end,1:end-1]

    a = alpha / Gamma
    c = -inv(Gamma)
    m = a[end,:]
  end
return RosenbrockStruct(
  nStage,
  alpha,
  Gamma,
  b,
  a,
  c,
  gamma,
  m,
  )
end
