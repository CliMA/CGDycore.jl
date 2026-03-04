import CGDycore: Integration as INT, IMEXRosenbrock as IR

FT = Float64
RK = INT.RungeKuttaMethod{FT}("SSP(5,3)")
ROS = INT.RosenbrockMethod{FT}("SSP-Knoth")

  # Y1 = y
  # k1 = h*f(Y1) + h*J(gamma11*k1)
  # Y2 = y + alpha21*k1
  # k2 = h*f(Y2) + h*J(gamma21*k1+gamma22*k2)
  # Y3 = y + alpha31*k1 + alpha32*k2
  # k3 = h*f(Y3) + h*J(gamma31*k1+gamma32*k2+gamma33*k3)
  # Y4 = y + alpha41*k1 + alpha42*k2 + alpha43*k3
  # k4 = h*f(Y4) + h*J(gamma41*k1+gamma42*k2+gamma43*k3+gamma44*k4)
  # yNew = b1*k1 + b2*k2 + b3*k3 + b4 * k4

  nStage = ROS.nStage + 1

  alphaR = [ROS.alpha[2:end,1:end]
            ROS.b']

  AE = [zeros(FT,1,nStage)
        alphaR zeros(FT,nStage-1)]
  bE = AE[end,:]

  AI = [zeros(FT,1,nStage)
        zeros(FT,nStage-1) alphaR * ROS.gamma * inv(alphaR)]
  bI = AI[end,:]
 
IMEX = IR.RosenbrockToIMEXDirk(ROS) 

ROSR = IR.IMEXDirkToRosenbrock(IMEX)

IMEXBo = IR.IMEXDirkMethod{FT}("Boscarino") 
ROSBo = IR.IMEXDirkToRosenbrock(IMEXBo)
OrderROSBo = IR.OrderConditionsRosenbrockW(ROSBo,3)
IMEXBHR = IR.IMEXDirkMethod{FT}("BHR(5,5,3)")
ROSBHR = IR.IMEXDirkToRosenbrock(IMEXBHR)
OrderROSBHR = IR.OrderConditionsRosenbrockW(ROSBHR,3)
@show OrderROSBHR

IMEXARS222 = IR.IMEXDirkMethod{FT}("ARS222") 
ROSARS222 = IR.IMEXDirkToRosenbrock(IMEXARS222)
OrderROSARS222 = IR.OrderConditionsRosenbrockW(ROSARS222,2)

O, ROSOptim = IR.FindRosenbrockMethod()

RKBo = INT.RungeKuttaMethod{FT}(ROSBo)
IR.StabilityRegion(RKBo)

RKBHR = INT.RungeKuttaMethod{FT}(ROSBHR)
IR.StabilityRegion(RKBHR)

RKROSOptim = INT.RungeKuttaMethod{FT}(ROSOptim)
IR.StabilityRegion(RKROSOptim)

