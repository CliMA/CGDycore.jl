function RungeKuttaExplicit!(time,V,dt,Fcn,CG,Global,Param)
  RK=Global.RK;
  f=Global.Cache.f
  Vn=Global.Cache.Vn

  @. Vn = V
  @inbounds for iStage=1:RK.nStage
    V .= Vn;
    @inbounds for jStage=1:iStage-1
      @views @. V = V + dt * RK.ARKE[iStage,jStage] * f[:,:,:,jStage]
    end
    Fcn(view(f,:,:,:,iStage),V,time + RK.cRKE[iStage] * dt,CG,Global,Param);
  end
  @. V = Vn;
  @inbounds for iStage=1:RK.nStage
    @views @. V = V + dt * RK.bRKE[iStage] * f[:,:,:,iStage]
  end
end

function RungeKuttaExplicit!(V,dt,Fcn,CG,Global,Param)
  RK=Global.RK;
  f=Global.Cache.f
  Vn=Global.Cache.Vn

  @. Vn = V
  @inbounds for iStage=1:RK.nStage
    V .= Vn;
    @inbounds for jStage=1:iStage-1
      @views @. V = V + dt * RK.ARKE[iStage,jStage] * f[:,:,:,jStage]
    end
#   Fcn(view(f,:,:,:,iStage),V,time + RK.cRKE[iStage] * dt,CG,Global,Param);
    @views Fcn(f[:,:,:,iStage],V,CG,Global,Param);
  end
  @. V = Vn;
  @inbounds for iStage=1:RK.nStage
    @views @. V = V + dt * RK.bRKE[iStage] * f[:,:,:,iStage]
  end
end

function RungeKuttaExplicitLS!(time,V,dt,Fcn,CG,Global,Param)
  RK=Global.RK;
  @views f=Global.Cache.f[:,:,:,1]
  Vn=Global.Cache.Vn

  @. Vn = V
  @inbounds for iStage=1:RK.nStage-1
    Global.Model.Param.time = time + RK.cRKE[iStage] * dt
    @time Fcn(f,V,CG,Global,Param);
    @. V = Vn + dt * RK.ARKE[iStage+1,iStage] * f
  end
  Fcn(f,V,CG,Global);
  @. V = Vn + dt * f
end

