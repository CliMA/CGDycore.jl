import CGDycore: Integration as INT, IMEXRosenbrock as IR

FT = Float64
IMEXTr = INT.IMEXDirkMethod{FT}("ARS343")
IMEX = IR.IMEXDirkMethod{FT}("ARS343")
AAE,AAI = IR.IMEXDirkNewtonFormulation(IMEX)


