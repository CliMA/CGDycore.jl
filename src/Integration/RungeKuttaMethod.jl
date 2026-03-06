function RungeKuttaExMethod{FT}() where FT<:AbstractFloat
  nStage = 0
  JacComp = false
  a = zeros(FT,0,0)
  b = zeros(FT,0)
  name = ""
  return RungeKuttaExMethod{FT}(
    name,
    nStage,
    a,
    b,
    JacComp
  )
end

function RungeKuttaExMethod{FT}(Ros::RosenbrockMethod) where FT<:AbstractFloat
  return RungeKuttaExMethod{FT}(
    "RK"*Ros.name,
    Ros.nStage,
    Ros.alpha,
    Ros.b,
    false
  )
end


function RungeKuttaExMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
  JacComp = false
  if str == "SSP(4,3)"
    nStage = 4
    a = [0 0 0 0
         1/2 0 0 0
         1/2 1/2 0 0
         1/6 1/6 1/6 0]
    b=[1/6;1/6;1/6;1/2]     
  elseif str == "SSP(5,3)"
    nStage = 5
    a= [0 0 0 0 0
        0.37726891511710 0 0 0 0
        0.37726891511710 0.37726891511710 0 0 0
        0.16352294089771 0.16352294089771 0.16352294089771 0 0
        0.14904059394856 0.14831273384724 0.14831273384724 0.34217696850008 0]
     b = [0.19707596384481; 0.11780316509765; 0.11709725193772; 0.27015874934251; 0.29786487010104]
  elseif str == "RK4"
    nStage = 4
    a = zeros(FT,nStage,nStage)
    a[2,1] = 1/2
    a[3,2] = 1/2
    a[4,3] = 1
    b = [1/6,1/3,1/3,1/6]
  end  

  return RungeKuttaExMethod{FT}(
    str,
    nStage,
    a,
    b,
    JacComp,
  )
end
