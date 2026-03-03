function FindRosenbrockMethod()
  FT = Float64
  RK = INT.RungeKuttaMethod{FT}("SSP(5,3)")
  Order = 3
  gammaD = 1.0
  Residual = RKRosenbrock()(RK,Order,gammaD)
  initial_x = ones(10)
  @views @. initial_x[1:3] = 0.0

# result = optimize(Residual, initial_x, BFGS())
  result = optimize(Residual, initial_x, NelderMead())
# result = optimize(Residual, zeros(10), SimulatedAnnealing())
  gamma = Optim.minimizer(result)

  Ros = INT.RosenbrockMethod{FT}(RK,gammaD,gamma)
  O = OrderConditionsRosenbrockW(Ros,Order)
  @show O
  return O,Ros
end  

