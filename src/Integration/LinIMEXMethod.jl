mutable struct LinIMEXStruct
  nStage::Int
  AHat::Array{Float64, 2}
  BHat::Array{Float64, 1}
  A::Array{Float64, 2}
  B::Array{Float64, 1}
end
function LinIMEXMethod()
  nStage=0
  AHat=zeros(0,0)
  BHat=zeros(0)
  A=zeros(0,0)
  B=zeros(0)
  return LinIMEXStruct(
    nStage,
    AHat,
    BHat,
    A,
    B,
  )
end

function LinIMEXMethod(Method)
str = Method
  if str == "ARS343"
    nStage = 4
    gamma = 0.4358665215084590
    b1 = -3/2 * gamma^2 + 4 * gamma - 1/4
    b2 =  3/2 * gamma^2 - 5 * gamma + 5/4
    a = [0 0 0 0;
         0 gamma 0 0;
        0 (1-gamma)/2 gamma 0;
        0 b1 b2 gamma]
    b = [0 b1 b2 gamma]

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
    aHat = [
      0 0 0 0;
      gamma 0 0 0;
      dA31 dA32 0 0;
      dA41 dA42 dA43 0]
    bHat = [  0 b1 b2 gamma]
  elseif str == "AR2"
    nStage = 3
    s2 = sqrt(2)
    a32 = 1/6*(3+2*s2)
    aHat = [0     0    0
            2-s2 0 0
            1-a32 a32 0]
    bHat = [1/(2*s2) 1/(2*s2) 1-1/s2]        
    
    a=[0 0 0
       1-1/s2 1-1/s2 0
       1/(2*s2) 1/(2*s2) 1-1/s2]
    b = [1/(2*s2) 1/(2*s2) 1-1/s2]        
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
  return LinIMEXStruct(
    nStage,
    AHat,
    BHat,
    A,
    B,
  )
end
