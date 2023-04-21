function FcnTracer!(F,U,time,CG,Global,Param)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    NumV,
    NumTr) = Global.Model

  dtau = 0.5 * Global.TimeStepper.dtauStage
  HorLimit = Global.Model.HorLimit
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  J = Global.Metric.J
  X = Global.Metric.X
  Temp1 = Global.Cache.Temp1
  @views DivTr = Global.Cache.Temp1[:,:,NumV+1:NumV+NumTr]
  @views JJ = Global.Cache.Temp1[:,:,NumV+NumTr+1]
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
  JJ .= 0.0
  x = StrideArray{Float64}(undef, 3)
  @show size(x)

  if HorLimit
    @views Limit!(qMin,qMax,U[:,:,NumV+1:NumV+NumTr],U[:,:,RhoPos],CG,Global)
  end  

  # Hyperdiffusion 
  @inbounds for iF in Global.Grid.BoundaryFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          JJ[iz,ind] += J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]
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
      @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTr[iz,ind,iT] += 0.5 * DivCG[iP,jP,iz] / CG.M[iz,ind]
          end
        end
      end
    end
  end

  @views ExchangeData3DSend(Temp1[:,:,1:NumV+NumTr+1],Global.Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          JJ[iz,ind] += J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]
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
      @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTr[iz,ind,iT] += 0.5 * DivCG[iP,jP,iz] / CG.M[iz,ind]
          end
        end
      end
    end
  end

  @views ExchangeData3DRecv!(Temp1[:,:,1:NumV+NumTr+1],Global.Exchange)

  @inbounds for iF in Global.Grid.BoundaryFaces
    if Param.StreamFun
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x1 = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x2 = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x3 = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            PsiLoc = fPsi(x1,x2,x3,time,Global,Param)
            PsiCG[iP,jP,iz] = PsiLoc
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

    @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

#   Tracer transport
    @inbounds for iT = 1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTrCG[iP,jP,iz] = DivTr[iz,ind,iT] / JJ[iz,ind]
          end
        end
      end
#     Hyperdiffusion, second Laplacian
      @views DivRhoGrad!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
        Global.Model.HyperDDiv)
      if Global.Model.Upwind
        @views DivUpwindRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],RhoCG,CG,
          Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],
          Global.ThreadCache)
      else
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
              Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)  
      end
    end


    
    if HorLimit
      @inbounds for iz=1:nz
        @inbounds for iT = 1:NumTr
          @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:]])
          @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:]])
          @views HorLimiter!(FCG[:,:,iz,iT+NumV],TrCG[:,:,iz,iT],RhoCG,RhoCG,dtau,
            Global.Metric.JC[:,:,iz,iF],CG.w,qMinS,qMaxS)
#         @views QP!(FCG[:,:,iz,iT+NumV],TrCG[:,:,iz,iT],RhoCG,RhoCG,dtau,
#           Global.Metric.JC[:,:,iz,iF],CG.w,qMinS,qMaxS)
        end
      end  
    end  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / JJ[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / JJ[iz,ind]
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
            x1 = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x2 = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x3 = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            PsiLoc = fPsi(x1,x2,x3,time,Global,Param)
            PsiCG[iP,jP,iz] = PsiLoc
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

    @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

#   Tracer transport
    @inbounds for iT = 1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTrCG[iP,jP,iz] = DivTr[iz,ind,iT] / JJ[iz,ind]
          end
        end
      end
#     Hyperdiffusion, second Laplacian
      @views DivRhoGrad!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
        Global.Model.HyperDDiv)
      if Global.Model.Upwind
        @views DivUpwindRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],RhoCG,CG,
          Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],
          Global.ThreadCache)
      else
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
              Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)  
      end
    end

    if HorLimit
      @inbounds for iz=1:nz
        @inbounds for iT = 1:NumTr
          @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:]])
          @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:]])
          @views HorLimiter!(FCG[:,:,iz,iT+NumV],TrCG[:,:,iz,iT],RhoCG,RhoCG,dtau,
            Global.Metric.JC[:,:,iz,iF],CG.w,qMinS,qMaxS)
#         @views QP!(FCG[:,:,iz,iT+NumV],TrCG[:,:,iz,iT],RhoCG,RhoCG,dtau,
#           Global.Metric.JC[:,:,iz,iF],CG.w,qMinS,qMaxS)
        end
      end
    end  
      
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / JJ[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / JJ[iz,ind]
          end
        end
      end  
    end
  end  

  ExchangeData3DRecv!(F,Global.Exchange)

end

function FcnTracerConv!(F,U,time::Float64,CG,Global,Param)
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
  x = StrideArray{Float64}(undef, StaticInt(3))
  x[1] = 1.0
  
  @inbounds for iF in Global.Grid.BoundaryFaces
    if Param.StreamFun
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x[1] = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x[2] = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x[3] = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            PsiLoc = @gc_preserve fPsi(x,time,Global,Param)
            PsiCG[iP,jP,iz] = PsiLoc
          end
        end
      end
      @views Curl!(v1CG,v2CG,PsiCG,CG,Global.Metric.dXdxI[:,:,:,:,:,:,iF],
        Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @. wCG = 0.0
    else
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x[1] = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x[2] = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x[3] = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            uu,vv = @gc_preserve fVel(x,time,Global,Param)
            v1CG[iP,jP,iz] = uu
            v2CG[iP,jP,iz] = vv
          end
        end
      end
      @. wCG = 0.0
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
    @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

    @inbounds for iT=1:NumTr
      @views DivConvRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
              Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += 0.5 * FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += 0.5 * FCG[iP,jP,iz,iT+NumV] / CG.MMass[iz,ind]
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
            x[1] = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x[2] = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x[3] = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            PsiLoc = @gc_preserve fPsi(x,time,Global,Param)
            PsiCG[iP,jP,iz] = PsiLoc
          end
        end
      end
      @views Curl!(v1CG,v2CG,PsiCG,CG,Global.Metric.dXdxI[:,:,:,:,:,:,iF],
        Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    else
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x[1] = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x[2] = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x[3] = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            uu,vv = @gc_preserve fVel(x,time,Global,Param)
            uu = 1.0
            vv = 1.0
            v1CG[iP,jP,iz] = uu
            v2CG[iP,jP,iz] = vv
          end
        end
      end
      @. wCG = 0.0
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
    @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

    @inbounds for iT = 1 : NumTr
      @views DivConvRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += 0.5 * FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += 0.5 * FCG[iP,jP,iz,iT+NumV] / CG.MMass[iz,ind]
          end
        end
      end
    end
  end

  ExchangeData3DRecv!(F,Global.Exchange)
end
