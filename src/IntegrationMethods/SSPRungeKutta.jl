function SSPRungeKutta!(time,V,dt,Fcn,CG,Metric,Phys,Cache,Global,Param,Profile)
  SSP=Global.TimeStepper.SSP
  nStage=SSP.nStage
  fV=Cache.fV
  fS=Cache.fS
  fRhoS=Cache.fRhoS
  VS=Cache.VS
  RhoS=Cache.RhoS
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr
  RhoPos=Global.Model.RhoPos

  @views @. RhoS[:,:,1] = V[:,:,RhoPos]
  if NumTr>0
    @views @. VS[:,:,:,1] = V[:,:,NumV+1:end]
  end  
  @inbounds for iStage = 1:nStage
    Global.TimeStepper.dtauStage = dt * SSP.ms[iStage]
    Fcn(fV,V,time + SSP.c[iStage] * dt,CG,Metric,Phys,Cache,Global,Param,Profile)
    @views @. fRhoS[:,:,iStage] = fV[:,:,RhoPos]
    if NumTr > 0
      @views @. fS[:,:,:,iStage] = fV[:,:,NumV+1:end]
    end
    @views @. RhoS[:,:,iStage+1] = 0.0 
    @views @. VS[:,:,:,iStage+1] = 0.0 
    @inbounds for jStage = 1:iStage
      if SSP.beta[iStage,jStage] > 0
        @views @. RhoS[:,:,iStage+1] += (SSP.alpha[iStage,jStage] * RhoS[:,:,jStage] +
          dt * SSP.beta[iStage,jStage] * fRhoS[:,:,jStage])
      else
        @views @. RhoS[:,:,iStage+1] += SSP.alpha[iStage,jStage] * RhoS[:,:,jStage]
      end    
      if NumTr>0
        if SSP.beta[iStage,jStage] > 0
          @views @. VS[:,:,:,iStage+1] += (SSP.alpha[iStage,jStage] * VS[:,:,:,jStage] +
            dt * SSP.beta[iStage,jStage] * fS[:,:,:,jStage])
        else
          @views @. VS[:,:,:,iStage+1] += SSP.alpha[iStage,jStage] * VS[:,:,:,jStage]
        end    
      end
    end  
    @views @. V[:,:,RhoPos] = RhoS[:,:,iStage+1]
    @views @. V[:,:,NumV+1:end] = VS[:,:,:,iStage+1]
  end
end

function SSPRungeKuttaGPU!(time,V,dt,Fcn,CG,Metric,Cache,Global,Param)
  SSP=Global.TimeStepper.SSP
  nStage=SSP.nStage
  fV=Cache.fV
  fS=Cache.fS
  fRhoS=Cache.fRhoS
  VS=Cache.VS
  RhoS=Cache.RhoS
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr
  RhoPos=Global.Model.RhoPos

  @views @. RhoS[:,:,1] = V[:,:,RhoPos]
  if NumTr>0
    @views @. VS[:,:,:,1] = V[:,:,NumV+1:end]
  end  
  @inbounds for iStage = 1:nStage
    Global.TimeStepper.dtauStage = dt * SSP.ms[iStage]
    Fcn(fV,V,time + SSP.c[iStage] * dt,CG,Metric,Cache,Global,Param)
    @views @. fRhoS[:,:,iStage] = fV[:,:,RhoPos]
    if NumTr > 0
      @views @. fS[:,:,:,iStage] = fV[:,:,NumV+1:end]
    end
    @views @. RhoS[:,:,iStage+1] = 0.0 
    @views @. VS[:,:,:,iStage+1] = 0.0 
    @inbounds for jStage = 1:iStage
      if SSP.beta[iStage,jStage] > 0
        @views @. RhoS[:,:,iStage+1] += (SSP.alpha[iStage,jStage] * RhoS[:,:,jStage] +
          dt * SSP.beta[iStage,jStage] * fRhoS[:,:,jStage])
      else
        @views @. RhoS[:,:,iStage+1] += SSP.alpha[iStage,jStage] * RhoS[:,:,jStage]
      end    
      if NumTr>0
        if SSP.beta[iStage,jStage] > 0
          @views @. VS[:,:,:,iStage+1] += (SSP.alpha[iStage,jStage] * VS[:,:,:,jStage] +
            dt * SSP.beta[iStage,jStage] * fS[:,:,:,jStage])
        else
          @views @. VS[:,:,:,iStage+1] += SSP.alpha[iStage,jStage] * VS[:,:,:,jStage]
        end    
      end
    end  
    @views @. V[:,:,RhoPos] = RhoS[:,:,iStage+1]
    @views @. V[:,:,NumV+1:end] = VS[:,:,:,iStage+1]
  end
end

