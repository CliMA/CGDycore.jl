mutable struct RosenbrockStruct{FT<:AbstractFloat}
  nStage::Int
  a::Array{FT, 2}
  c::Array{FT, 2}
  gamma::FT
  m::Array{FT, 1}
end

function RosenbrockStruct{FT}() where FT<:AbstractFloat
  nStage = 0
  a = zeros(FT,0,0)
  c = zeros(FT,0,0)
  gamma = FT(0)
  m = zeros(FT,0)
  return RosenbrockStruct{FT}(
    nStage,
    a,
    c,
    gamma,
    m,
  )
end

function RosenbrockStruct{FT}(Method) where FT<:AbstractFloat
  str = Method
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
    aCPU = alpha / Gamma
    cCPU = -inv(Gamma)
    mCPU = Gamma'\b
    aCPU=[aCPU
       mCPU']

  elseif str == "RosEul"
    nStage = 1
    alpha = zeros(FT,nStage,nStage)
    b = zeros(FT,nStage)
    b[1] = 1
    Gamma = zeros(FT,nStage,nStage)
    Gamma[1,1] = 1
    gamma = FT(1)        
    aCPU = alpha / Gamma
    cCPU = -inv(Gamma)
    mCPU = Gamma'\b
    aCPU=[aCPU
       mCPU']
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
    aCPU = alpha / Gamma
    cCPU = -inv(Gamma)
    mCPU = Gamma'\b
    aCPU=[aCPU
       mCPU']
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

  elseif str == "ROS34PRW"
    nStage = 4
    alpha = zeros(nStage,nStage)
    alpha[2,1] = 8.7173304301691801e-01       
    alpha[3,1] = 1.4722022879435914e+00       
    alpha[3,2] = -3.1840250568090289e-01      
    alpha[4,1] = 8.1505192016694938e-01       
    alpha[4,2] = 5.0000000000000000e-01       
    alpha[4,3] = -3.1505192016694938e-01      

    Gamma = zeros(nStage,nStage)
    Gamma[1,1] = 4.3586652150845900e-01 
    Gamma[2,1] = -8.7173304301691801e-01
    Gamma[2,2] = 4.3586652150845900e-01 
    Gamma[3,1] = -1.2855347382089872e+00
    Gamma[3,2] = 5.0507005541550687e-01
    Gamma[3,3] =  4.3586652150845900e-01
    Gamma[4,1] = -4.8201449182864348e-01
    Gamma[4,2] = 2.1793326075422950e-01
    Gamma[4,3] = -1.7178529043404503e-01
    Gamma[4,4] =  4.3586652150845900e-01
    gamma = 4.3586652150845900e-01

    b = zeros(nStage)
    b[1] = 3.3303742833830591e-01
    b[2] = 7.1793326075422947e-01  
    b[3] = -4.8683721060099439e-01
    b[4] = 4.3586652150845900e-01

    aCPU = alpha / Gamma
    cCPU = -inv(Gamma)
    mCPU = Gamma'\b
    aCPU=[aCPU
       mCPU']
  #=
    ROS34PRW.
    γ = 4.3586652150845900e-01
    α21 = 8.7173304301691801e-01 γ21 = -8.7173304301691801e-01
    α31 = 1.4722022879435914e+00 γ31 = -1.2855347382089872e+00
    α32 = -3.1840250568090289e-01 γ32 = 5.0507005541550687e-01
    α41 = 8.1505192016694938e-01 γ41 = -4.8201449182864348e-01
    α42 = 5.0000000000000000e-01 γ42 = 2.1793326075422950e-01
    α43 = -3.1505192016694938e-01 γ43 = -1.7178529043404503e-01
    b1 = 3.3303742833830591e-01
    b2 = 7.1793326075422947e-01
    b3 = -4.8683721060099439e 01
    b4 = 4.3586652150845900e-01
  =#
  elseif str == "ROSWRODAS3"
    nStage = 4
    alpha = zeros(nStage,nStage)
    alpha = [0    0     0   0
             0    0     0   0
             1.   0     0   0
             0.75 -0.25 0.5 0]
    Gamma = [0.5     0       0       0 
             1.      0.5     0       0 
             -0.25   -0.25   0.5     0
             1. / 12 1. / 12 -2. / 3 0.5]
    b = [5. / 6, -1. / 6, -1. / 6, 0.5]
    gamma = 0.5
    aCPU = alpha / Gamma
    cCPU = -inv(Gamma)
    mCPU = Gamma'\b
    aCPU=[aCPU
       mCPU']
  elseif str == "ROS3P"
#=
    nStage = 3
    aCPU = zeros(nStage,nStage)
    aCPU[2,1] = 1.26794919243112
    aCPU[3,1] = 1.26794919243112
    cCPU = zeros(nStage,nStage)
    cCPU[2,1] = -1.607695154586736
    cCPU[3,1] = -3.464101615137755
    cCPU[3,2] = -1.732050807568877
    gamma = 
=#

  end
  return RosenbrockStruct{FT}(
    nStage,
    aCPU,
    cCPU,
    gamma,
    mCPU,
  )
end
