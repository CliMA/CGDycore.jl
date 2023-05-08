function Fcn!(F,U,CG,Global,Param,::Val{:TestGrad})

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model

  Phys=Global.Phys    
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  J = Global.Metric.J
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views JJ = Global.Cache.Temp1[:,:,NumV+NumTr+1]
  @views JRho = Global.Cache.Temp1[:,:,NumV+NumTr+2]
  @views JRhoF = Global.Cache.Temp1[:,:,NumV+NumTr+3]
  FCG=Global.Cache.FCC
  FwCG=Global.Cache.FwCC
  RhoCG = Global.Cache.RhoCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  zPG = Global.Cache.zPG

  F .= 0.0
  JJ .= 0.0
  JRho .= 0.0
  JRhoF .= 0.0

  @inbounds for iF = 1 : Global.Grid.NumFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          JJ[iz,ind] += J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]
          JRho[iz,ind] += (J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]) 
        end
        @inbounds for iz=1:nz-1
          JRhoF[iz,ind] += (J[iP,jP,2,iz,iF] + J[iP,jP,1,iz+1,iF])
        end
      end
    end
  end

  @inbounds for iF = 1 : Global.Grid.NumFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind]
        end
      end
    end
    @. FCG = 0.0
    @. FwCG = 0.0
    @views Gradient!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos]
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FwCG[iP,jP,iz+1] 
        end
      end  
    end
  end  

  @views @. F[:,:,uPos] /= JJ
  @views @. F[:,:,vPos] /= JJ
  @views @. F[1:nz-1,:,wPos] /= JRhoF[1:nz-1,:]

end

function Fcn!(F,U,CG,Global,Param,::Val{:TestFunCGrad})

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model

  Phys=Global.Phys    
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  J = Global.Metric.J
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views JJ = Global.Cache.Temp1[:,:,NumV+NumTr+1]
  @views JRho = Global.Cache.Temp1[:,:,NumV+NumTr+2]
  @views JRhoF = Global.Cache.Temp1[:,:,NumV+NumTr+3]
  FCG=Global.Cache.FCC
  FwCG=Global.Cache.FwCC
  cCG = Global.Cache.RhoCG
  fCG = Global.Cache.ThCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  zPG = Global.Cache.zPG

  F .= 0.0
  JJ .= 0.0
  JRho .= 0.0
  JRhoF .= 0.0

  @inbounds for iF = 1 : Global.Grid.NumFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          JJ[iz,ind] += J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]
          JRho[iz,ind] += (J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]) 
        end
        @inbounds for iz=1:nz-1
          JRhoF[iz,ind] += (J[iP,jP,2,iz,iF] + J[iP,jP,1,iz+1,iF])
        end
      end
    end
  end

  @inbounds for iF = 1 : Global.Grid.NumFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          cCG[iP,jP,iz] = U[iz,ind,1]
          fCG[iP,jP,iz] = U[iz,ind,2]
        end
      end
    end
    @. FCG = 0.0
    @. FwCG = 0.0
    @views FunCGradient!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],cCG,fCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos]
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FwCG[iP,jP,iz+1] 
        end
      end  
    end
  end  

  @views @. F[:,:,uPos] /= JJ
  @views @. F[:,:,vPos] /= JJ
  @views @. F[1:nz-1,:,wPos] /= JRhoF[1:nz-1,:]

end

