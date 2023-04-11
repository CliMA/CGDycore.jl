function FcnTracer!(F,U,time,CG,Global,Param)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    NumV,
    NumTr) = Global.Model

  dtau = Global.TimeStepper.dtauStage
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  X = Global.Metric.X
  Temp1 = Global.Cache.Temp1
  @views DivTr = Global.Cache.Temp1[:,:,NumV+1:NumV+NumTr]
  DivCG=Global.Cache.DivC
  FCG=Global.Cache.FCC
  DivTrCG=Global.Cache.DivC
  RhoCG = Global.Cache.RhoCG
  ThCG = Global.Cache.ThCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  PsiCG = Global.Cache.wCCG
  TrCG = Global.Cache.TrCG
  DivTr .= 0.0
  F .= 0.0
  qMin = Global.Cache.qMin
  qMax = Global.Cache.qMax

  @views Limit!(qMin,qMax,U[:,:,NumV+1:NumV+NumTr],U[:,:,RhoPos],CG,Global)

  # Hyperdiffusion 
  @inbounds for iF in Global.Grid.BoundaryFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
        end
      end
    end
    @inbounds for iT=1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            ThCG[iP,jP,iz] = U[iz,ind,iT+NumV]
          end
        end
      end
      @views FDivGrad2VecDSS!(DivTrCG,ThCG,RhoCG,CG,Global,iF)
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTr[iz,ind,iT] += DivCG[iP,jP,iz] / CG.M[iz,ind]
          end
        end
      end
    end
  end

  @views ExchangeData3DSend(Temp1[:,:,1:NumV+NumTr],Global.Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
        end
      end
    end
    @inbounds for iT=1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            ThCG[iP,jP,iz] = U[iz,ind,iT+NumV]
          end
        end
      end
      @views FDivGrad2VecDSS!(DivTrCG,ThCG,RhoCG,CG,Global,iF)
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTr[iz,ind,iT] += DivCG[iP,jP,iz] / CG.M[iz,ind]
          end
        end
      end
    end
  end

  ExchangeData3DRecv!(Temp1[:,:,1:NumV+NumTr],Global.Exchange)

  @inbounds for iF in Global.Grid.BoundaryFaces
    if Param.StreamFun
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x = SVector{3}(0.5 * (X[iP,jP,1,:,iz,iF] .+ X[iP,jP,2,:,iz,iF]))
            @views PsiCG[iP,jP,iz] = fPsi(x,time,Global,Param)
          end
        end  
      end  
      @views Curl!(v1CG,v2CG,PsiCG,CG,Global.Metric.dXdxI[:,:,:,:,:,:,iF],
        Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    else
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x = SVector{3}(0.5 * (X[iP,jP,1,:,iz,iF] .+ X[iP,jP,2,:,iz,iF]))
            @views (uu,vv) = fVel(x,time,Global,Param)
            v1CG[iP,jP,iz] = uu
            v2CG[iP,jP,iz] = vv
            wCG[iP,jP,iz+1] = 0.0
          end
        end  
      end  
    end
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,iT+NumV]
          end  
        end
      end
    end
    @. FCG = 0.0

    @views BoundaryW!(wCG[:,:,:],v1CG[:,:,:],v2CG[:,:,:],CG,
       Global.Metric.J,Global.Metric.dXdxI[:,:,:,1,:,:,iF])
    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF)

#   Tracer transport
    @inbounds for iT = 1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTrCG[iP,jP,iz] = DivTr[iz,ind,iT]
          end
        end
      end
#     Hyperdiffusion, second Laplacian
      @views FDivRhoGrad2Vec!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG,Global,iF)
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
              Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)  
      end
    end

    @inbounds for iz=1:nz
      @inbounds for iT = 1:NumTr
        @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:]])
        @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:]])
        @views SecantQ!(FCG[:,:,iz,iT+NumV],TrCG[:,:,iz,iT],RhoCG,RhoCG,dtau,
          Global.Metric.JC[:,:,iz,iF],CG.w,qMinS,qMaxS)
      end
    end  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
      end  
    end
  end  
  
  ExchangeData3DSend(F,Global.Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    if Param.StreamFun
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x = SVector{3}(0.5 * (X[iP,jP,1,:,iz,iF] .+ X[iP,jP,2,:,iz,iF]))
            @views PsiCG[iP,jP,iz] = fPsi(x,time,Global,Param)
          end 
        end  
      end  
      @views Curl!(v1CG,v2CG,PsiCG,CG,Global.Metric.dXdxI[:,:,:,:,:,:,iF],
        Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    else
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x = SVector{3}(0.5 * (X[iP,jP,1,:,iz,iF] .+ X[iP,jP,2,:,iz,iF]))
            @views (uu,vv) = fVel(x,time,Global,Param)
            v1CG[iP,jP,iz] = uu
            v2CG[iP,jP,iz] = vv
            wCG[iP,jP,iz+1] = 0.0 
          end 
        end  
      end  
    end 
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,iT+NumV]
          end  
        end
      end
    end
    @. FCG = 0.0

#   BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views BoundaryW!(wCG[:,:,:],v1CG[:,:,:],v2CG[:,:,:],CG,
       Global.Metric.J,Global.Metric.dXdxI[:,:,:,1,:,:,iF])
    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF)

#   Tracer transport
    @inbounds for iT = 1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTrCG[iP,jP,iz] = DivTr[iz,ind,iT]
          end
        end
      end
#     Hyperdiffusion, second Laplacian
      @views FDivRhoGrad2Vec!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG,Global,iF)
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
              Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)  
      end
    end

    @inbounds for iz=1:nz
      @inbounds for iT = 1:NumTr
        @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:]])
        @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:]])
        @views SecantQ!(FCG[:,:,iz,iT+NumV],TrCG[:,:,iz,iT],RhoCG,RhoCG,dtau,
          Global.Metric.JC[:,:,iz,iF],CG.w,qMinS,qMaxS)
      end
    end  
      
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
      end  
    end
  end  

  ExchangeData3DRecv!(F,Global.Exchange)

end

