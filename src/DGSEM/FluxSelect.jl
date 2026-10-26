################################################################################
# Selection of the volume (two-point) fluxes and the interface flux for the
# split-form DGSEM (FcnSplit!). Sets Model.FluxAverageH, Model.FluxAverageV,
# Model.RiemannSolver and Model.NonConservativeFlux.
################################################################################

function SetDensityMean(DensityMeanDG)
  if DensityMeanDG == "Arith"
    return ArithMean()
  elseif DensityMeanDG == "Log"
    return LogMean()
  elseif DensityMeanDG == "Gamma"
    return GammaMean()
  else
    error("DensityMeanDG = $DensityMeanDG, choose Arith, Log or Gamma")
  end
end

function SetFluxes!(Model,DG,FT,FluxDG,DensityMeanDG,InterfaceFluxDG,Param,Phys,pAuxPos,GPAuxPos,cS)
  RhoPos = Model.RhoPos
  uPos = Model.uPos
  vPos = Model.vPos
  wPos = Model.wPos
  RhoThPos = Model.RhoThPos

  # The volume fluxes return (rho, rho u, rho v, rho w, rho theta) in this order
  if (RhoPos,uPos,vPos,wPos,RhoThPos) != (1,2,3,4,5) || Model.NumV != 5
    error("split-form volume fluxes need the variables (Rho,Rhou,Rhov,Rhow,RhoTh) at positions 1:5")
  end

  if FluxDG == "KennedyGruberGrav"
    FluxAverage = KennedyGruberGrav()(Val(pAuxPos),Val(GPAuxPos))
  elseif FluxDG == "PotTempTEC"
    FluxAverage = PotTempTEC()(Val(pAuxPos),Val(GPAuxPos),SetDensityMean(DensityMeanDG),FT(Phys.Gamma))
  elseif FluxDG == "PotTempEC"
    FluxAverage = PotTempEC()(Val(pAuxPos),Val(GPAuxPos),SetDensityMean(DensityMeanDG),FT(Phys.Gamma))
  elseif FluxDG == "PotTempETEC"
    FluxAverage = PotTempETEC()(Val(pAuxPos),Val(GPAuxPos),FT(Phys.Gamma))
  else
    error("FluxDG = $FluxDG is not available for the split-form DGSEM")
  end
  Model.FluxAverageH = FluxAverage
  Model.FluxAverageV = FluxAverage

  if InterfaceFluxDG == "RiemannLMARS"
    Model.RiemannSolver = RiemannLMARS()(Val(RhoPos),Val(uPos),Val(vPos),Val(wPos),
      Val(RhoThPos),Val(1),FT(cS))
  elseif InterfaceFluxDG == "RiemannExLMARS"
    Model.RiemannSolver = RiemannExLMARS()(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,1)
  elseif InterfaceFluxDG == "RiemannExnerLMARS"
    Model.RiemannSolver = RiemannExnerLMARS()(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,1)
  elseif InterfaceFluxDG == "RiemannExPLMARS"
    Model.RiemannSolver = RiemannExPLMARS()(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,1)
  else
    error("InterfaceFluxDG = $InterfaceFluxDG is not available for the split-form DGSEM")
  end

  Model.NonConservativeFlux = BuoyancyFlux()(RhoPos,GPAuxPos)
end
