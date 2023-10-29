mutable struct SSPRungeKuttaStruct
  nStage::Int
  alpha::Array{Float64, 2}
  beta::Array{Float64, 2}
  c::Array{Float64, 1}
  ms::Array{Float64, 1}
end  
function SSPRungeKuttaMethod()
 nStage=0
 alpha=zeros(0,0)
 beta=zeros(0,0)
 c=zeros(0)
 ms=zeros(0)
 return SSPRungeKuttaStruct(
   nStage,
   alpha,
   beta,
   c,
   ms,
   )
end
function SSPRungeKuttaMethod(alpha::Array{Float64, 2},beta::Array{Float64, 2})
 c=zeros(size(alpha,1)) 
 for i = 2:size(alpha,1)
   c[i] = sum(beta[i-1,:])  
 end  
 ms=zeros(0)
 return SSPRungeKuttaStruct(
   size(alpha,1),
   alpha,
   beta,
   c,
   ms,
   )
end 

function SSPRungeKuttaMethod(Method)
  if Method == "SSP-Knoth"
    nStage = 3  
    alpha = [1 0 0
             3/4 1/4 0
             1/3 0 2/3]
    beta = [1 0 0
            0 1/4 0
            0 0 2/3]
  elseif Method == "SSP32"
    nStage = 3  
    alpha = [1 0 0
             3/4 1/4 0
             1/3 0 2/3]
    beta = [1 0 0
            0 1/4 0
            0 0 2/3]
    c = [0,1,1/2]
  elseif Method == "SSPRK1"
    nStage = 1  
    alpha = zeros(1,1)
    alpha[1,1] = 1.0
    beta = zeros(1,1)
    beta[1,1] = 1.0
    c = [0]
  end          
  ms = zeros(nStage)
  for i = 1 : nStage
    for j = i : nStage
      if alpha[j,i] != 0.0
        ms[i] = max(ms[i], beta[j,i] / alpha[j,i])
      end
    end
  end  
 return SSPRungeKuttaStruct(
   nStage,
   alpha,
   beta,
   c,
   ms,
   )
end 
  

