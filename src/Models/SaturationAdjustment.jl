function SaturationAdjustment!(Rho,RhoTheta,RhoV,RhoC,Phys)

#
# Find dRhoV such that with
# RhoThetaA = RhoTheta + Lv * RhoTheta / T * dRhoV
# RhoVA = RhoV + dRhoV
# RhoCA = RhoC - dRhoV
# 
# Equilibrium(Rho,RhoThetaA,RhoVA,RhoCA) = 0


# Compute T
 T = fTemp(Rho,RhoV,RhoC,RhoTheta,Phys)
# Compute iinitial internal energy
 Lv = LatHeat(T,Phys)
 e = InternalEnergy(Rho,RhoV,RhoC,T,Phys)
# Newton-Method

# Compute RhoVS
  for iTer = 1 : 3
    pVS = fpvs1(T,Phys)
    RhoVS = pVS / (Phys.Rv *  T) 
    RhoVNew = min(RhoVS,RhoV + RhoC)
    RhoCNew = RhoV + RhoC - RhoVNew

    eNew = InternalEnergy(Rho,RhoVNew,RhoCNew,T,Phys)

    f = eNew - e
    dfdT = dInternalEnergydT(Rho,RhoVNew,RhoCNew,T,Phys) 
    @show f,dfdT
    TNew = T - f/dfdT
    @show T,TNew,f
    @show RhoV,RhoVNew
    @show RhoC,RhoCNew
    T = TNew
    RhoV = RhoVNew
    RhoC = RhoCNew
  end
  RhoTheta = fThetaV(Rho,RhoV,RhoC,T,Phys)
# Compute RhoThetaNew
end
