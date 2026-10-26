################################################################################
# Selection of the volume (two-point) fluxes and the interface flux for the
# split-form DGSEM (FcnSplit!). Sets Model.FluxAverageH, Model.FluxAverageV,
# Model.RiemannSolver and Model.NonConservativeFlux.
################################################################################

function SetFluxes!(Model,DG,FT,FluxDG,InterfaceFluxDG,Param,Phys,pAuxPos,GPAuxPos,cS)
  RhoPos = Model.RhoPos
  uPos = Model.uPos
  vPos = Model.vPos
  wPos = Model.wPos
  RhoThPos = Model.RhoThPos

  if FluxDG == "KennedyGruber"
    Model.FluxAverageH = KennedyGruber()(RhoPos,uPos,vPos,wPos,RhoThPos,pAuxPos,DG.H)
    Model.FluxAverageV = KennedyGruber()(RhoPos,uPos,vPos,wPos,RhoThPos,pAuxPos,DG.V)
  elseif FluxDG == "KennedyGruberGrav"
    Model.FluxAverageH = KennedyGruberGrav()(Val(RhoPos),Val(uPos),Val(vPos),Val(wPos),
      Val(RhoThPos),Val(pAuxPos),Val(GPAuxPos),DG.H)
    Model.FluxAverageV = KennedyGruberGrav()(Val(RhoPos),Val(uPos),Val(vPos),Val(wPos),
      Val(RhoThPos),Val(pAuxPos),Val(GPAuxPos),DG.V)
  else
    error("FluxDG = $FluxDG is not available for the split-form DGSEM")
  end

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
