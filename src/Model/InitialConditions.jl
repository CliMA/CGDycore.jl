function InitialConditions(CG,Global,Param)  
  Model = Global.Model
  nz = Global.Grid.nz
  NumV = Model.NumV
  NumTr = Model.NumTr
  U = zeros(Float64,nz,CG.NumG,NumV+NumTr)
  U[:,:,Model.RhoPos]=Project(fRho,0.0,CG,Global,Param)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=ProjectVec(fVel,0.0,CG,Global,Param)
  if Global.Model.Thermo == "InternalEnergy"
    U[:,:,Model.ThPos]=Project(fIntEn,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  elseif Global.Model.Thermo == "TotalEnergy"
    U[:,:,Model.ThPos]=Project(fTotEn,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  else
    U[:,:,Model.ThPos]=Project(fTheta,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end
  if NumTr>0
    U[:,:,Model.RhoVPos+Model.NumV]=Project(fQv,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end
  Global.pBGrd = Project(fpBGrd,0.0,CG,Global,Param)
  Global.RhoBGrd = Project(fRhoBGrd,0.0,CG,Global,Param)
  return U
end  

function InitialConditionsAdvection(CG,Global,Param)
  Model = Global.Model
  nz = Global.Grid.nz
  NumV = Model.NumV
  NumTr = Model.NumTr
  U = zeros(Float64,nz,CG.NumG,NumV+NumTr)
  U[:,:,Model.RhoPos]=Project(fRho,0.0,CG,Global,Param)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=ProjectVec(fVel,0.0,CG,Global,Param)
  U[:,:,Model.wPos]=ProjectW(fVelW,0.0,CG,Global,Param)
  for i = 1 : NumTr
    Model.ProfTr="AdvectionSphereDCMIPQ1"
    @views U[:,:,Model.NumV+i]=Project(fTr,0.0,CG,Global,Param).*U[:,:,Model.RhoPos]
  end  
  return U
end  