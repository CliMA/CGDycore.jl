function RungeKuttaMethod{FT}() where FT<:AbstractFloat
  nStage = 0
  a = zeros(FT,0,0)
  b = zeros(FT,0)
  name = ""
  return RungeKuttaMethod{FT}(
    name,
    nStage,
    a,
    b,
  )
end

function RungeKuttaMethod{FT}(Ros::RosenbrockMethod) where FT<:AbstractFloat
  return RungeKuttaMethod{FT}(
    "RK"*Ros.name,
    Ros.nStage,
    Ros.alpha,
    Ros.b,
  )
end


function RungeKuttaMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
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
  end  

  return RungeKuttaMethod{FT}(
    str,
    nStage,
    a,
    b,
  )
end
