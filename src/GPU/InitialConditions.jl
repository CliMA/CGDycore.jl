function InitialConditions(backend,FTB,CG,Metric,Phys,Global,Profile,Param)
  Model = Global.Model
  Nz = Global.Grid.nz
  NF = Global.Grid.NumFaces
  NumV = Model.NumV
  NumTr = Model.NumTr
  N = CG.OrdPoly + 1
  Glob = CG.Glob
  X = Metric.X
  time = 0


  # Ranges
  NzG = min(div(256,N*N),Nz)
  group = (N * N, NzG, 1)
  ndrange = (N * N, Nz, NF)

  U = KernelAbstractions.zeros(backend,FTB,Nz,CG.NumG,NumV+NumTr)
  @views Rho = U[:,:,Model.RhoPos]
  @views u = U[:,:,Model.uPos]
  @views v = U[:,:,Model.vPos]
  @views w = U[:,:,Model.wPos]
  @views RhoTh = U[:,:,Model.ThPos]
  KRhoFunCKernel! = RhoFunCKernel!(backend, group)
  KRhoThFunCKernel! = RhoThFunCKernel!(backend, group)
  KuvwFunCKernel! = uvwFunCKernel!(backend, group)

  KRhoFunCKernel!(Profile,Rho,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  KuvwFunCKernel!(Profile,u,v,w,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  KRhoThFunCKernel!(Profile,RhoTh,time,Glob,X,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  if Model.RhoVPos > 0
    @views RhoV = U[:,:,Model.RhoVPos]
    KRhoVFunCKernel! = RhoVFunCKernel!(backend, group)
    KRhoVFunCKernel!(Profile,RhoV,time,Glob,X,ndrange=ndrange)
    KernelAbstractions.synchronize(backend)
  end  
  if Model.RhoCPos > 0
    @views RhoC = U[:,:,Model.RhoCPos]
    KRhoCFunCKernel! = RhoCFunCKernel!(backend, group)
    KRhoCFunCKernel!(Profile,RhoC,time,Glob,X,ndrange=ndrange)
    KernelAbstractions.synchronize(backend)
  end  

#=
  if Global.Model.Profile
    Profile = TestRes(Global.Phys)
  else
    Profile = zeros(0)  
  end  
  U = zeros(Float64,nz,CG.NumG,NumV+NumTr)
  U[:,:,Model.RhoPos]=Project(fRho,0.0,CG,Metric,Phys,Global,Param,Profile)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=ProjectVec(fVel,0.0,CG,Metric,Phys,Global,Param)
  if Global.Model.Thermo == "InternalEnergy"
    U[:,:,Model.ThPos]=Project(fIntEn,0.0,CG,Metric,Phys,Global,Param,Profile).*U[:,:,Model.RhoPos]
  elseif Global.Model.Thermo == "TotalEnergy"
    U[:,:,Model.ThPos]=Project(fTotEn,0.0,CG,Metric,Phys,Global,Param,Profile).*U[:,:,Model.RhoPos]
  else
    U[:,:,Model.ThPos]=Project(fTheta,0.0,CG,Metric,Phys,Global,Param,Profile).*U[:,:,Model.RhoPos]
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
      U[:,:,Model.RhoVPos+Model.NumV]=Project(fQv,0.0,CG,Metric,Phys,Global,Param,Profile).*U[:,:,Model.RhoPos]
    end
    if Model.RhoCPos > 0  
      U[:,:,Model.RhoCPos+Model.NumV]=Project(fQc,0.0,CG,Metric,Phys,Global,Param,Profile).*U[:,:,Model.RhoPos]
    end
  end
  if Global.Model.ModelType == "Conservative"
    @views @. U[:,:,Model.uPos] *= U[:,:,Model.RhoPos]  
    @views @. U[:,:,Model.vPos] *= U[:,:,Model.RhoPos]  
  end  
  Global.pBGrd = Project(fpBGrd,0.0,CG,Metric,Phys,Global,Param,Profile)
  Global.RhoBGrd = Project(fRhoBGrd,0.0,CG,Metric,Phys,Global,Param,Profile)
  Global.ThetaBGrd = zeros(nz,CG.NumG)
  Global.TBGrd = zeros(nz,CG.NumG)
  if Global.Model.Geos
    (Global.UGeo,Global.VGeo) = ProjectVec(fVelGeo,0.0,CG,Metric,Phys,Global,Param)
  end  
  =#
  return U
end  

function InitialConditionsAdvection(backend,FTB,CG,Metric,Phys,Global,Profile,Param)
  Model = Global.Model
  Nz = Global.Grid.nz
  NF = Global.Grid.NumFaces
  NumV = Model.NumV
  NumTr = Model.NumTr
  N = CG.OrdPoly + 1
  Glob = CG.Glob
  X = Metric.X
  time = 0

  # Ranges
  NzG = min(div(256,N*N),Nz)
  group = (N * N, NzG, 1)
  ndrange = (N * N, Nz, NF)

  @show NF

  U = KernelAbstractions.zeros(backend,FTB,Nz,CG.NumG,NumV+NumTr)
  @views Rho = U[:,:,Model.RhoPos]
  @views u = U[:,:,Model.uPos]
  @views v = U[:,:,Model.vPos]
  @views w = U[:,:,Model.wPos]
  @views Tr = U[:,:,Model.NumV+1]
  KRhoFunCKernel! = RhoFunCKernel!(backend, group)
  KRhoFunCKernel!(Profile,Rho,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  KuvwFunCKernel! = uvwFunCKernel!(backend, group)
  KuvwFunCKernel!(Profile,u,v,w,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  KTrFunCKernel! = TrFunCKernel!(backend, group)
  KTrFunCKernel!(Profile,Tr,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  return U
end  
