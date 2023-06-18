function InitialConditions(CG,Global,Param)  
  Model = Global.Model
  nz = Global.Grid.nz
  NumV = Model.NumV
  NumTr = Model.NumTr
  NumG = CG.NumG
  OrdPoly = CG.OrdPoly

  if Global.Model.Profile
    Profile = TestRes(Global.Phys)
  else
    Profile = zeros(0)  
  end  
  U = zeros(Float64,nz,CG.NumG,NumV+NumTr)
  U[:,:,Model.RhoPos]=Project(fRho,0.0,CG,Global,Param,Profile)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=ProjectVec(fVel,0.0,CG,Global,Param)
  if Global.Model.Thermo == "InternalEnergy"
    U[:,:,Model.ThPos]=Project(fIntEn,0.0,CG,Global,Param,Profile).*U[:,:,Model.RhoPos]
  elseif Global.Model.Thermo == "TotalEnergy"
    U[:,:,Model.ThPos]=Project(fTotEn,0.0,CG,Global,Param,Profile).*U[:,:,Model.RhoPos]
  else
    U[:,:,Model.ThPos]=Project(fTheta,0.0,CG,Global,Param,Profile).*U[:,:,Model.RhoPos]
  end
  if Model.PertTh
    for iF = 1 : size(U[:,:,Model.ThPos],2)    
      @. U[:,iF,Model.ThPos] *= (1.0 + 1.e-5 * (2.0*rand() - 1.0) * (400.0 <= Global.Metric.zP[:,iF] <= 2000.0))    
    end
    J = Global.Metric.J
    ThLoc = zeros(nz,NumG)
    for iF=1:Global.Grid.NumFaces
      for iz=1:nz
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            ind = CG.Glob[i,j,iF]
            ThLoc[iz,ind] += U[iz,ind,Model.ThPos] * (J[i,j,1,iz,iF] + J[i,j,2,iz,iF]) / CG.M[iz,ind]
          end
        end
      end
    end
    ExchangeData!(ThLoc,Global.Exchange)
    @. U[:,:,Model.ThPos] = ThLoc
  end  
  if NumTr>0
    if Model.RhoVPos > 0  
      U[:,:,Model.RhoVPos+Model.NumV]=Project(fQv,0.0,CG,Global,Param,Profile).*U[:,:,Model.RhoPos]
    end
    if Model.RhoCPos > 0  
      U[:,:,Model.RhoCPos+Model.NumV]=Project(fQc,0.0,CG,Global,Param,Profile).*U[:,:,Model.RhoPos]
    end
  end
  if Global.Model.ModelType == "Conservative"
    @views @. U[:,:,Model.uPos] *= U[:,:,Model.RhoPos]  
    @views @. U[:,:,Model.vPos] *= U[:,:,Model.RhoPos]  
  end  
  Global.pBGrd = Project(fpBGrd,0.0,CG,Global,Param,Profile)
  Global.RhoBGrd = Project(fRhoBGrd,0.0,CG,Global,Param,Profile)
  Global.ThetaBGrd = zeros(nz,CG.NumG)
  Global.TBGrd = zeros(nz,CG.NumG)
  if Global.Model.Geos
    (Global.UGeo,Global.VGeo) = ProjectVec(fVelGeo,0.0,CG,Global,Param)
  end  
  return U
end  

function InitialConditionsAdvection(CG,Global,Param)
  Model = Global.Model
  nz = Global.Grid.nz
  NumV = Model.NumV
  NumTr = Model.NumTr
  U = zeros(Float64,nz,CG.NumG,NumV+NumTr)
  Profile = zeros(0)
  U[:,:,Model.RhoPos]=Project(fRho,0.0,CG,Global,Param,Profile)
  for i = 1 : NumTr
    @views U[:,:,Model.NumV+i]=Project(fTr,0.0,CG,Global,Param,Profile).*U[:,:,Model.RhoPos]
  end  
  return U
end  
