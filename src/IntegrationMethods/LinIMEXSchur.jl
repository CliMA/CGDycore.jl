function LinIMEXSchur!(V,dt,Fcn,Jac,CG,Global,Param)
  LinIMEX=Global.LinIMEX;
  nV1=size(V,1);
  nV2=size(V,2);
  nV3=size(V,3);
  nJ=nV1*nV2*nV3;
  nStage=LinIMEX.nStage;
  f=Global.Cache.f
  fV=Global.Cache.fV
  Vn=Global.Cache.Vn
  Ymyn=Global.Cache.Ymyn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr

  J = Global.J
  @. Vn = V
  @inbounds for iStage = 2 : nStage
    @views Fcn(f[:,:,:,iStage-1],V,CG,Global,Param);
    if iStage == 2
      Jac(J,V,CG,Global,Param)
    end  
    @views @. fV = LinIMEX.AHat[iStage,iStage-1] * f[:,:,:,iStage-1]
    @inbounds for jStage = 1 : iStage - 2
      if LinIMEX.AHat[iStage,jStage] != 0.0
        @views @. fV += LinIMEX.AHat[iStage,jStage]*f[:,:,:,jStage] 
      end
    end
    @inbounds for jStage = 1 : iStage - 2
      if LinIMEX.A[iStage,jStage+1] != 0.0
        @views @. fV -= (LinIMEX.A[iStage,jStage+1] / dt) * Ymyn[:,:,:,jStage] 
      end
    end
    J.CompTri=true
    @views SchurSolve!(Ymyn[:,:,:,iStage-1],fV,J,dt/LinIMEX.A[iStage,iStage],Global);
    @views @. V = Ymyn[:,:,:,iStage-1] + Vn
  end
  if LinIMEX.BHat[nStage] != 0.0
    @views Fcn(f[:,:,:,nStage],V,CG,Global,Param);
    @views @. V = Vn + dt * LinIMEX.BHat[nStage] * f[:,:,:,nStage]
  else
    @. V = Vn  
  end  
  @inbounds for iStage = 1 : nStage - 1
    if LinIMEX.BHat[iStage] != 0.0
      @views @. V += (dt * LinIMEX.BHat[iStage]) * f[:,:,:,iStage] 
    end
  end
  @inbounds for iStage = 1 : nStage - 1
    if LinIMEX.B[iStage+1] != 0.0
      @views @. V += LinIMEX.B[iStage+1] * Ymyn[:,:,:,iStage]
    end
  end
end


function LinIMEXSchur!(V,dt,Fcn,FcnE,FcnI,Jac,CG,Global,Param)
  LinIMEX=Global.LinIMEX;
  nV1=size(V,1);
  nV2=size(V,2);
  nV3=size(V,3);
  nJ=nV1*nV2*nV3;
  nStage=LinIMEX.nStage;
  f=Global.Cache.f
  fV=Global.Cache.fV
  Vn=Global.Cache.Vn
  Ymyn=Global.Cache.Ymyn
  NumV=Global.Model.NumV
  NumTr=Global.Model.NumTr

  J = Global.J
  @. Vn = V
  fE = similar(V)
  fI = similar(V)
  @inbounds for iStage = 2 : nStage
    @views Fcn(f[:,:,:,iStage-1],V,CG,Global,Param);
    FcnE(fE,V,CG,Global,Param);
    FcnI(fI,V,CG,Global,Param);
    @views @. f[:,:,:,iStage-1] = fE + fI
    if iStage == 2
      Jac(J,V,CG,Global,Param)
    end  
    @views @. fV = LinIMEX.AHat[iStage,iStage-1] * f[:,:,:,iStage-1]
    @inbounds for jStage = 1 : iStage - 2
      if LinIMEX.AHat[iStage,jStage] != 0.0
        @views @. fV += LinIMEX.AHat[iStage,jStage]*f[:,:,:,jStage] 
      end
    end
    @inbounds for jStage = 1 : iStage - 2
      if LinIMEX.A[iStage,jStage+1] != 0.0
        @views @. fV -= (LinIMEX.A[iStage,jStage+1] / dt) * Ymyn[:,:,:,jStage] 
      end
    end
    J.CompTri=true
    @views SchurSolve!(Ymyn[:,:,:,iStage-1],fV,J,dt/LinIMEX.A[iStage,iStage],Global);
    @views @. V = Ymyn[:,:,:,iStage-1] + Vn
  end
  if LinIMEX.BHat[nStage] != 0.0
    @views Fcn(f[:,:,:,nStage],V,CG,Global,Param);
    @views @. V = Vn + dt * LinIMEX.BHat[nStage] * f[:,:,:,nStage]
  else
    @. V = Vn  
  end  
  @inbounds for iStage = 1 : nStage - 1
    if LinIMEX.BHat[iStage] != 0.0
      @views @. V += (dt * LinIMEX.BHat[iStage]) * f[:,:,:,iStage] 
    end
  end
  @inbounds for iStage = 1 : nStage - 1
    if LinIMEX.B[iStage+1] != 0.0
      @views @. V += LinIMEX.B[iStage+1] * Ymyn[:,:,:,iStage]
    end
  end
end


