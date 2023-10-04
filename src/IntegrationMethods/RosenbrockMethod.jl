mutable struct RosenbrockStruct{FT<:AbstractFloat,
                               AT1<:AbstractArray,
                               AT2<:AbstractArray}
  nStage::Int
  alpha::AT2
  Gamma::AT2
  b::AT1
  a::AT2
  c::AT2
  gamma::Float64
  m::AT1
end

function RosenbrockStruct{FT}(backend) where FT<:AbstractFloat
  nStage = 0
  alpha = KernelAbstractions.zeros(backend,FT,0,0)
  Gamma = KernelAbstractions.zeros(backend,FT,0,0)
  b = KernelAbstractions.zeros(backend,FT,0)
  a = KernelAbstractions.zeros(backend,FT,0,0)
  c = KernelAbstractions.zeros(backend,FT,0,0)
  gamma = 0
  m = KernelAbstractions.zeros(backend,FT,0)
  d = KernelAbstractions.zeros(backend,FT,0)
  return RosenbrockStruct{FT,
                          typeof(b),
                          typeof(alpha)}(
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

function RosenbrockStruct{FT}(backend,Method) where FT<:AbstractFloat
  str = Method
  if str == "SSP-Knoth"
    nStage = 3
    alpha = KernelAbstractions.zeros(backend,FT,nStage,nStage)
    alpha[2,1] = 1
    alpha[3,1] = 1/4
    alpha[3,2] = 1/4
    b = KernelAbstractions.zeros(backend,FT,nStage)
    b[1] = 1/6
    b[2] = 1/6
    b[3] = 2/3
    Gamma = KernelAbstractions.zeros(backend,FT,nStage,nStage)
    Gamma[1,1] = 1
    Gamma[2,2] = 1
    Gamma[3,1] = -3/4
    Gamma[3,2] = -3/4
    Gamma[3,3] = 1
    gamma = 1        
    a = alpha / Gamma
    c = -inv(Gamma)
    m = Gamma'\b
    a=[a
       m']

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
return RosenbrockStruct{FT,
                        typeof(b),
                        typeof(alpha)}(
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
