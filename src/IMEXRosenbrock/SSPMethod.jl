mutable struct SSPMethod{FT<:AbstractFloat}
  nStage::Int
  alpha::Array{FT, 2}
  beta::Array{FT, 2}
  a::Array{FT, 2}
  b::Array{FT, 1}
end

function SSPMethod{FT}() where FT<:AbstractFloat
  nStage = 0
  alpha = zeros(FT,0,0)
  beta = zeros(FT,0,0)
  a = zeros(FT,0,0)
  b = zeros(FT,0)
  return SSPMethod{FT}(
    nStage,
    alpha,
    beta,
    a,
    b,
  )
end


function SSPMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
  if str == "SSP(4,3)"
    nStage = 4
    alpha = zeros(FT,nStage+1,nStage+1)
    beta = zeros(FT,nStage+1,nStage+1)
    b = zeros(FT,nStage)


    a = [0 0 0 0
         1/2 0 0 0
         1/2 1/2 0 0
         1/6 1/6 1/6]

    b=[1/6;1/6;1/6;1/2]     

    gamma = zeros(FT,nStage,nStage)
    gamma[1,1] = 1
    gamma[2,2] = 1
    gamma[3,1] = -3/4
    gamma[3,2] = -3/4
    gamma[3,3] = 1
    gammaD = FT(1)        

  elseif str == "6S4O"
#   6S4O.C/W-method ? Using LDDRK coefficients given in [17] with D 0:25:
    nStage = 6
    a21 = 0.28878526699679
    a31 = 0.10893125722541
    a32 = 0.27283594644263
    a41 = 0.10893125722541
    a42 = 0.13201701492152
    a43 = 0.47167254854945
    a51 = 0.10893125722541
    a52 = 0.13201701492152
    a53 = 0.38911623225517
    a54 = 0.06600540453183
    a61 = 0.10893125722541
    a62 = 0.13201701492152
    a63 = 0.38911623225517
    a64 = -0.59203884581148
    a65 = 0.79248022128095

    g11 = 0.25
    g21 = -0.45345741148076
    g31 = -0.34182832909418
    g32 = 0.00000000000000
    g41 = -1.93637949137395
    g42 = 0.62221779527294
    g43 = 0.83345812222713
    g51 = -1.10275049376267
    g52 = 0.47337577919072
    g53 = 0.27833333985558
    g54 = -0.02663940566679
    g61 = -0.97465070482040
    g62 = 0.04287310605107
    g63 = 0.98104398325919
    g64 = 0.59370081382312
    g65 = -0.97639882505842

    b1 = 0.10893125722541
    b2 = 0.13201701492152
    b3 = 0.38911623225517
    b4 = -0.59203884581148
    b5 = 0.47385028714844
    b6 = 0.48812405426094
    alpha = [0   0   0   0   0   0
      a21 0   0   0   0   0
      a31 a32 0   0   0   0
      a41 a42 a43 0   0   0
      a51 a52 a53 a54 0   0
      a61 a62 a63 a64 a65 0
      b1  b2  b3  b4  b5  b6]
    gamma = [g11 0   0   0   0   0
      g21 g11 0   0   0   0
      g31 g32 g11 0   0   0
      g41 g42 g43 g11 0   0
      g51 g52 g53 g54 g11 0
      g61 g62 g63 g64 g65 g11]
    gammaD = g11  
    b = [b1; b2; b3; b4; b5; b6]  
  end  
  a = alpha / gamma
  c = -inv(gamma)
  m = gamma'\b
  return SSPMethod{FT}(
    nStage,
    alpha,
    beta,
    a,
    b,
  )
end
