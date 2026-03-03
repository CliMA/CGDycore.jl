abstract type Residual end
Base.@kwdef struct RKRosenbrock <: Residual end

function (residual::RKRosenbrock)(RK,Order,gammaD)
  function local_residual(x)
    FT = eltype(x)
    Ros = INT.RosenbrockMethod{FT}(RK,gammaD,x)
    O = OrderConditionsRosenbrockW(Ros,Order)
    f = norm(O)
    return f
  end
  return local_residual
end
