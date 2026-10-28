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

# Conservative interface fluxes are the same for both sides
function TwoSided(RiemannSolver!)
  @inline function RiemannTwoSided!(FL,FR,VLL,VRR,AuxL,AuxR,n1,n2,n3)
    RiemannSolver!(FL,VLL,VRR,AuxL,AuxR,n1,n2,n3)
    @unroll for iv = 1 : length(FL)
      FR[iv] = FL[iv]
    end
  end
  return RiemannTwoSided!
end

# Selection in three steps: discretization (DG, CG), formulation (ModelType
# Conservative, VectorInvariant), then volume and interface flux. FluxDG or
# InterfaceFluxDG = "Default" selects the default flux of the formulation.
const DefaultFluxes = Dict(
  "Conservative" => ("KennedyGruberGrav", "RiemannLMARS"),
  "VectorInvariant" => ("VectorInvariant", "RiemannVectorInvariant"))

function SetFluxes!(Model,Discretization,FT,FluxDG,DensityMeanDG,InterfaceFluxDG,Param,Phys,pAuxPos,GPAuxPos,cS)
  # CG has no two-point volume and interface fluxes
  Discretization == "DG" || return

  ModelType = Model.ModelType
  if !haskey(DefaultFluxes, ModelType)
    error("ModelType = $ModelType is not available for DG, choose $(join(keys(DefaultFluxes), " or "))")
  end
  FluxDGDef, InterfaceFluxDGDef = DefaultFluxes[ModelType]
  FluxDG = FluxDG == "Default" ? FluxDGDef : FluxDG
  InterfaceFluxDG = InterfaceFluxDG == "Default" ? InterfaceFluxDGDef : InterfaceFluxDG
  @info "DG fluxes" ModelType FluxDG InterfaceFluxDG

  if ModelType == "Conservative"
    SetFluxesConservative!(Model,FT,FluxDG,DensityMeanDG,InterfaceFluxDG,Param,Phys,pAuxPos,GPAuxPos,cS)
  elseif ModelType == "VectorInvariant"
    SetFluxesVectorInvariant!(Model,FT,FluxDG,InterfaceFluxDG,Phys,pAuxPos,GPAuxPos,cS)
  end
end

# Conservative split form (rho, rho v, rho theta)
function SetFluxesConservative!(Model,FT,FluxDG,DensityMeanDG,InterfaceFluxDG,Param,Phys,pAuxPos,GPAuxPos,cS)
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
    error("FluxDG = $FluxDG is not available for ModelType Conservative, choose KennedyGruberGrav, PotTempTEC, PotTempEC or PotTempETEC")
  end
  Model.FluxAverageH = FluxAverage
  Model.FluxAverageV = FluxAverage

  if InterfaceFluxDG == "RiemannLMARS"
    Model.RiemannSolver = TwoSided(RiemannLMARS()(Val(RhoPos),Val(uPos),Val(vPos),Val(wPos),
      Val(RhoThPos),Val(1),FT(cS)))
  elseif InterfaceFluxDG == "RiemannExLMARS"
    Model.RiemannSolver = TwoSided(RiemannExLMARS()(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,1))
  elseif InterfaceFluxDG == "RiemannExnerLMARS"
    Model.RiemannSolver = TwoSided(RiemannExnerLMARS()(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,1))
  elseif InterfaceFluxDG == "RiemannExPLMARS"
    Model.RiemannSolver = TwoSided(RiemannExPLMARS()(Param,Phys,RhoPos,uPos,vPos,wPos,RhoThPos,1))
  else
    error("InterfaceFluxDG = $InterfaceFluxDG is not available for ModelType Conservative, choose RiemannLMARS, RiemannExLMARS, RiemannExnerLMARS or RiemannExPLMARS")
  end

  Model.NonConservativeFlux = BuoyancyFlux()(RhoPos,GPAuxPos)
end

# Vector invariant form (rho, v, rho theta), arXiv:2610.05820. The Exner
# pressure is kept in Aux at ExnerAuxPos.
const ExnerAuxPos = 3

function SetFluxesVectorInvariant!(Model,FT,FluxDG,InterfaceFluxDG,Phys,pAuxPos,GPAuxPos,cS)
  if (Model.RhoPos,Model.uPos,Model.vPos,Model.wPos,Model.RhoThPos) != (1,2,3,4,5) || Model.NumV != 5
    error("VectorInvariant needs the variables (Rho,u,v,w,RhoTh) at positions 1:5")
  end
  if Model.NumAux < ExnerAuxPos || ExnerAuxPos in (pAuxPos,GPAuxPos)
    error("VectorInvariant keeps the Exner pressure at Aux position $ExnerAuxPos, use NumAux >= $ExnerAuxPos")
  end

  if Model.ThetaMeanDG == "Log"
    ThMean = LogMean()
  elseif Model.ThetaMeanDG == "Arith"
    ThMean = ArithMean()
  else
    error("ThetaMeanDG = $(Model.ThetaMeanDG), choose Log or Arith")
  end
  if Model.SoundSpeedDG == "Local"
    cFix = zero(FT)
  elseif Model.SoundSpeedDG == "Fixed"
    cFix = FT(cS)
  else
    error("SoundSpeedDG = $(Model.SoundSpeedDG), choose Local or Fixed")
  end

  if FluxDG == "VectorInvariant"
    FluxAverage = VectorInvariantFlux()(Val(ExnerAuxPos),Val(GPAuxPos),FT(Phys.Cpd),ThMean)
  else
    error("FluxDG = $FluxDG is not available for ModelType VectorInvariant, use VectorInvariant")
  end
  Model.FluxAverageH = FluxAverage
  Model.FluxAverageV = FluxAverage

  if InterfaceFluxDG == "RiemannVectorInvariant"
    Model.RiemannSolver = RiemannVectorInvariant()(Val(pAuxPos),Val(ExnerAuxPos),Val(GPAuxPos),
      FT(Phys.Cpd),FT(Phys.Gamma),ThMean,cFix)
  else
    error("InterfaceFluxDG = $InterfaceFluxDG is not available for ModelType VectorInvariant, use RiemannVectorInvariant")
  end
end
