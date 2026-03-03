mutable struct MISMethod{FT<:AbstractFloat} <: IntegrationMethod
  nStage::Int
  beta::Array{FT, 2}
  alpha::Array{FT, 2}
  gamma::Array{FT, 2}
  d::Array{FT, 1}
  FastMethod::IntegrationMethod
end


function MISMethod{FT}() where FT<:AbstractFloat
  nStage = 0
  beta = zeros(FT,0,0)
  alpha = zeros(FT,0,0)
  gamma = zeros(FT,0,0)
  d = zeros(FT, 0)
  FastMethod = NoMethod()
  return MISMethod{FT}(
    nStage,
    beta,
    alpha,
    gamma,
    d,
    FastMethod,
  )
end

function MISMethod{FT}(Method) where FT<:AbstractFloat
  FastMethod = NoMethod()
  str = Method
  if str == "MISRK4"
    nStage = 4
	    beta = zeros(FT, nStage+1, nStage)
            alpha = zeros(FT, nStage+1, nStage)
            gamma = zeros(FT, nStage+1, nStage)
            d = zeros(FT, nStage+1, 1)
    	    beta[2, 1] = 0.38758444641450318
            beta[3, 1] = -2.5318448354142823E-002
            beta[3, 2] = 0.38668943087310403
            beta[4, 1] = 0.20899983523553325
            beta[4, 2] = -0.45856648476371231
            beta[4, 3] = 0.43423187573425748
            beta[5, 1] = -0.10048822195663100
            beta[5, 2] = -0.46186171956333327
            beta[5, 3] = 0.83045062122462809
            beta[5, 4] = 0.27014914900250392
    
            alpha[3, 2] = 0.52349249922385610
            alpha[4, 2] = 1.1683374366893629
            alpha[4, 3] = -0.75762080241712637
            alpha[5, 2] = -3.6477233846797109E-002
            alpha[5, 3] = 0.56936148730740477
            alpha[5, 4] = 0.47746263002599681
    
            gamma[3, 2] = 0.13145089796226542
            gamma[4, 2] = -0.36855857648747881
            gamma[4, 3] = 0.33159232636600550
            gamma[5, 2] = -6.5767130537473045E-002
            gamma[5, 3] = 4.0591093109036858E-002
            gamma[5, 4] = 6.4902111640806712E-002
    
		d[2] = beta[2, 1]
		d[3] = beta[3, 1] + beta[3, 2]
		d[4] = beta[4, 1] + beta[4, 2] + beta[4, 3]
		d[5] = beta[5, 1] + beta[5, 2] + beta[5, 3] + beta[5, 4]
  elseif str == "RK3"
    nStage = 3
    beta = zeros(FT, nStage+1, nStage)
    alpha = zeros(FT, nStage+1, nStage)
    gamma = zeros(FT, nStage+1, nStage)
    d = zeros(FT, nStage+1, 1)
    beta[2, 1] = 1/3
    beta[3, 2] = 1/2
    beta[4, 3] = 1
    d[2] = 1/3
    d[3] = 1/2
    d[4] = 1
  elseif str == "RKJeb"
    nStage = 3
    beta = zeros(FT, nStage+1, nStage)
    alpha = zeros(FT, nStage+1, nStage)
    gamma = zeros(FT, nStage+1, nStage)
    d = zeros(FT, nStage+1)

    beta[2,1] = 2.0492941060709863e-001
    beta[3,1] = -4.5477553356788974e-001
    beta[3,2] = 9.5613538239378981e-001
    beta[4,1] = -3.5970281266252929e-002
    beta[4,2] = -1.5363649484946584e-001
    beta[4,3] = 7.0259062712330234e-001

    gamma[3,2] = -8.2176071248067006e-001
    gamma[4,2] = -3.8080670922635063e-001
    gamma[4,3] = 4.5653105107801978e-001

    alpha[3,2] = 7.0302371060435331e-001
    alpha[4,2] = 4.2492220536139252e-001
    alpha[4,3] = 5.4545718243573982e-001

    d[2] = beta[2, 1]
    d[3] = beta[3, 1] + beta[3, 2]
    d[4] = beta[4, 1] + beta[4, 2] + beta[4, 3]
   end
  return MISMethod{FT}(
    nStage,
    beta,
    alpha,
    gamma,
    d,
    FastMethod,
    )
end

mutable struct MISSemiMethod{FT<:AbstractFloat} <: IntegrationMethod
  nStage::Int
  beta::Array{FT, 2}
  alpha::Array{FT, 2}
  gamma::Array{FT, 2}
  d::Array{FT, 1}
  FastMethod::IntegrationMethod
end
function MISSemiMethod{FT}() where FT<:AbstractFloat
  MIS = MISMethod{FT}()
  return MISSemiMethod{FT}(
    MIS.nStage,
    MIS.beta,
    MIS.alpha,
    MIS.gamma,
    MIS.d,
    MIS.FastMethod,
    )
end  
function MISSemiMethod{FT}(Method) where FT<:AbstractFloat
  MIS = MISMethod{FT}(Method)
  return MISSemiMethod{FT}(
    MIS.nStage,
    MIS.beta,
    MIS.alpha,
    MIS.gamma,
    MIS.d,
    MIS.FastMethod,
    )
end  

mutable struct MISLinMethod{FT<:AbstractFloat} <: IntegrationMethod
  nStage::Int
  beta::Array{FT, 2}
  alpha::Array{FT, 2}
  gamma::Array{FT, 2}
  d::Array{FT, 1}
  FastMethod::IntegrationMethod
end
function MISLinMethod{FT}() where FT<:AbstractFloat
  MIS = MISMethod{FT}()
  return MISLinMethod{FT}(
    MIS.nStage,
    MIS.beta,
    MIS.alpha,
    MIS.gamma,
    MIS.d,
    MIS.FastMethod,
    )
end  
function MISLinMethod{FT}(Method) where FT<:AbstractFloat
  MIS = MISMethod{FT}(Method)
  return MISLinMethod{FT}(
    MIS.nStage,
    MIS.beta,
    MIS.alpha,
    MIS.gamma,
    MIS.d,
    MIS.FastMethod,
    )
end  

