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

  elseif str == "ROSRK3"
    nStage = 3
    alpha = KernelAbstractions.zeros(backend,FT,nStage,nStage)
    alpha[2,1] = 1/3
    alpha[3,2] = 1/2
    b = KernelAbstractions.zeros(backend,FT,nStage)
    b[3] = 1
    gamma = 1
    Gamma = KernelAbstractions.zeros(backend,FT,nStage,nStage)
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
  end
return RosenbrockStruct{FT}(
  nStage,
  aCPU,
  cCPU,
  gamma,
  mCPU,
  )
end
