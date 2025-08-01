NVTX.@annotate function Rosenbrock!(V,dt,Fcn!,Jac,CG,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  ROS = Global.TimeStepper.ROS
  nStage = ROS.nStage
  k = Cache.k
  fV = Cache.fV
  Vn = Cache.Vn
  Global.TimeStepper.dtauStage = dt  # Oswald

  JCache.CompTri = true
  @. Vn = V
  @inbounds for iStage = 1 : nStage
    @. V = Vn
    @inbounds for jStage = 1 : iStage-1
      @views @. V = V + ROS.a[iStage,jStage] * k[:,:,:,jStage]
    end
    Fcn!(fV,V,CG,Metric,Phys,Cache,Exchange,Global,Param,DiscType)
    if iStage == 1
      Jac(JCache,V,CG,Metric,Phys,Cache,Global,Param,DiscType)
    end  
    @inbounds for jStage = 1 : iStage - 1
      fac = ROS.c[iStage,jStage] / dt
      @views @. fV = fV + fac * k[:,:,:,jStage]
    end
    @views Solve!(k[:,:,:,iStage],fV,JCache,dt*ROS.gamma,Cache,Global)
  end
  @. V = Vn
  @inbounds for iStage = 1 : nStage
    @views @. V = V + ROS.m[iStage] * k[:,:,:,iStage]
  end
end

function RosenbrockDG!(V,dt,Fcn!,Jac,CG,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,DiscType)

  ROS = Global.TimeStepper.ROS
  nStage = ROS.nStage
  k = Cache.k
  fV = Cache.fV
  Vn = Cache.Vn
  Global.TimeStepper.dtauStage = dt  # Oswald

  JCache.CompTri = true
  @. Vn = V
  @inbounds for iStage = 1 : nStage
    @. V = Vn
    @inbounds for jStage = 1 : iStage-1
      @views @. V = V + ROS.a[iStage,jStage] * k[:,:,:,:,jStage]
    end
    FcnPrepare!(V,CG,Metric,Phys,Cache,Exchange,Global,Param,DiscType)
    Fcn!(fV,V,CG,Metric,Phys,Cache,Exchange,Global,Param,DiscType)
    if iStage == 1
      Jac(JCache,V,CG,Metric,Phys,Cache,Global,Param,DiscType)
    end  
    @inbounds for jStage = 1 : iStage - 1
      fac = ROS.c[iStage,jStage] / dt
      @views @. fV = fV + fac * k[:,:,:,:,jStage]
    end
    @views Solve!(k[:,:,:,iStage],fV,JCache,dt*ROS.gamma,Cache,Global)
  end
  @. V = Vn
  @inbounds for iStage = 1 : nStage
    @views @. V = V + ROS.m[iStage] * k[:,:,:,:,iStage]
  end
end

function RosenbrockSchurMIS!(V,dt,Fcn,R,Jac,CG,Global,Param)
  ROS=Global.TimeStepper.ROS
  nV1=size(V,1)
  nV2=size(V,2)
  nV3=size(V,3)
  nJ=nV1*nV2*nV3
  nStage=ROS.nStage
  k=Cache.k
  fV=Cache.fV
  Vn=Cache.Vn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr

  J = Global.J
  J.CompTri=true
  @. Vn = V
  @inbounds for iStage=1:nStage
    @. V = Vn
    @inbounds for jStage=1:iStage-1
      @views @. V = V + ROS.a[iStage,jStage]*k[:,:,:,jStage]
    end
    Fcn(fV,V,CG,Global,Param)
    @. fV = fV + R
    if iStage == 1
      Jac(J,V,CG,Global,Param)
    end
    @inbounds for jStage=1:iStage-1
        @views @. fV = fV + (ROS.c[iStage,jStage]/dt)*k[:,:,:,jStage]
    end
    @views SchurSolve!(k[:,:,:,iStage],fV,J,dt*ROS.Gamma[iStage,iStage],Global)
  end
  @. V = Vn
  @inbounds for iStage=1:nStage
    @views @. V = V + ROS.m[iStage] * k[:,:,:,iStage]
  end

end

function RosenbrockDSchur!(V,dt,Fcn,Jac,CG,Global)
  ROS=Global.TimeStepper.ROS
  nV1=size(V,1)
  nV2=size(V,2)
  nV3=size(V,3)
  nJ=nV1*nV2*nV3
  nStage=ROS.nStage
  k=Cache.k
  fV=Cache.fV
  Vn=Cache.Vn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr

  J = Global.J
  Jac(J,V,CG,Global)
  @. Vn = V
  @inbounds for iStage=1:nStage
    @. V = Vn
    @inbounds for jStage=1:iStage-1
      @views @. V = V + ROS.a[iStage,jStage]*k[:,:,:,jStage]
    end
    Fcn(fV,V,CG,Global)
    @inbounds for jStage=1:iStage-1
        @views @. fV = fV + (ROS.c[iStage,jStage]/dt)*k[:,:,:,jStage]
    end
    J.CompTri=true
    @views SchurSolve!(k[:,:,:,iStage],fV,J,dt*ROS.Gamma[iStage,iStage],Global)
  end
  @. V = Vn
  @inbounds for iStage=1:nStage
    @views @. V[:,:,1:NumV+NumTr] = V[:,:,1:NumV+NumTr] + ROS.m[iStage] * k[:,:,:,iStage]
  end
end


function RosenbrockSSP!(V,dt,Fcn,Jac,CG,Global)
  ROS=Global.TimeStepper.ROS
  SSP=Global.TimeStepper.ROS.SSP
  nV1=size(V,1)
  nV2=size(V,2)
  nV3=size(V,3)
  nJ=nV1*nV2*nV3
  nStage=ROS.nStage
  k=Cache.k
  fV=Cache.fV
  fS=Cache.fS
  fRhoS=Cache.fRhoS
  VS=Cache.VS
  RhoS=Cache.RhoS
  Vn=Cache.Vn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr
  RhoPos=Global.Model.RhoPos

  J = Global.J
  J.CompTri=true
  @. Vn = V
  if NumTr>0
    @views @. VS[:,:,:,1] = V[:,:,NumV+1:end]
  end  
  @inbounds for iStage = 1:nStage
    Fcn(fV,V,CG,Global)
    if iStage == 1
      Jac(J,V,CG,Global)
    end  
    if NumTr > 0
      @views @. fS[:,:,:,iStage] = fV[:,:,NumV+1:end]
    end  
    @inbounds for jStage = 1:iStage-1
      @views @. fV[:,:,1:NumV] = fV[:,:,1:NumV] + (ROS.c[iStage,jStage]/dt)*k[:,:,:,jStage]
    end
    @views SchurSolve!(k[:,:,:,iStage],fV,J,dt*ROS.Gamma[iStage,iStage],Global)
    @views @. V[:,:,1:NumV] = Vn[:,:,1:NumV]
    @inbounds for jStage = 1:iStage
      @views @. V[:,:,1:NumV] = V[:,:,1:NumV] + ROS.a[iStage+1,jStage]*k[:,:,:,jStage]
    end
    if NumTr>0
      @views @. VS[:,:,:,iStage+1] = 0.0 
      @inbounds for jStage = 1:iStage
        if SSP.beta[iStage,jStage] > 0
          if Global.Model.HorLimit  
            Fac = SSP.beta[iStage,jStage]/SSP.alpha[iStage,jStage]
            @views HorLimiter!(VS[:,:,:,iStage+1],SSP.alpha[iStage,jStage],VS[:,:,:,jStage],Fac*dt,
              fS[:,:,:,jStage],V[:,:,RhoPos],CG,Global)            
          else 
            @views @. VS[:,:,:,iStage+1] += (SSP.alpha[iStage,jStage] * VS[:,:,:,jStage] +
              dt * SSP.beta[iStage,jStage] * fS[:,:,:,jStage])
          end  
        else
          @views @. VS[:,:,:,iStage+1] += SSP.alpha[iStage,jStage] * VS[:,:,:,jStage]
        end    
      end
      @views @. V[:,:,NumV+1:end] = VS[:,:,:,iStage+1]
    end  
  end
end

