function FcnTracer!(F,U,time,CG,Metric,Phys,Cache,Exchange,Global,Param,Profile)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    NumV,
    NumTr) = Global.Model

  @show minimum(U[:,:,NumV+1]),maximum(U[:,:,NumV+1])
  dtau = Global.TimeStepper.dtauStage
  HorLimit = Global.Model.HorLimit
  Grav=Phys.Grav    
  OP=CG.OrdPoly+1;
  DoF = CG.DoF
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  J = Metric.J
  X = Metric.X
  Temp1 = Cache.Temp1
  @views DivTr = Cache.Temp1[:,:,NumV+1:NumV+NumTr]
  DivCG=Cache.DivC
  FCG=Cache.FCC
  DivTrCG=Cache.DivC
  RhoCG = Cache.RhoCG
  ThCG = Cache.ThCG
  v1CG = Cache.v1CG
  v2CG = Cache.v2CG
  wCG = Cache.wCG
  PsiCG = Cache.wCCG
  TrCG = Cache.TrCG
  DivTr .= 0.0
  F .= 0.0
  qMin = Cache.qMin
  qMax = Cache.qMax

  if HorLimit
    @views Limit!(qMin,qMax,U[:,:,NumV+1:NumV+NumTr],U[:,:,RhoPos],CG,Global)
  end  

  # Hyperdiffusion 
  @inbounds for iF in Global.Grid.BoundaryFaces
    @inbounds for iD = 1 : DoF
      ind = CG.Glob[iD,iF]
      @inbounds for iz=1:nz
        RhoCG[iD,iz] = U[iz,ind,RhoPos]
      end
    end
    @inbounds for iT=1:NumTr
      @inbounds for iD = 1 : DoF
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          ThCG[iD,iz] = U[iz,ind,iT+NumV]
        end
      end
      @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
        Metric.dXdxI[:,:,:,:,:,iF],Metric.J[:,:,:,iF],Global.ThreadCache)
      @inbounds for iD = 1 : DoF
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          DivTr[iz,ind,iT] += DivCG[iD,iz] / CG.M[iz,ind]
        end
      end
    end
  end

  @views Parallels.ExchangeData3DSend(Temp1[:,:,1:NumV+NumTr],Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    @inbounds for iD = 1 : DoF
      ind = CG.Glob[iD,iF]
      @inbounds for iz=1:nz
        RhoCG[iD,iz] = U[iz,ind,RhoPos]
      end
    end
    @inbounds for iT=1:NumTr
      @inbounds for iD = 1 : DoF
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          ThCG[iD,iz] = U[iz,ind,iT+NumV]
        end
      end
      @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
        Metric.dXdxI[:,:,:,:,:,iF],Metric.J[:,:,:,iF],Global.ThreadCache)
      @inbounds for iD = 1 : DoF
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          DivTr[iz,ind,iT] += DivCG[iD,iz] / CG.M[iz,ind]
        end
      end
    end
  end

  @views Parallels.ExchangeData3DRecv!(Temp1[:,:,1:NumV+NumTr],Exchange)

  @inbounds for iF in Global.Grid.BoundaryFaces
    if Param.StreamFun
      @inbounds for iD = 1 : DoF
        @inbounds for iz=1:nz
          x1 = 0.5 * (X[iD,1,1,iz,iF] + X[iD,2,1,iz,iF])
          x2 = 0.5 * (X[iD,1,2,iz,iF] + X[iD,2,2,iz,iF])
          x3 = 0.5 * (X[iD,1,3,iz,iF] + X[iD,2,3,iz,iF])
          xS = SVector{3}(x1, x2 ,x3)
          PsiLoc = fPsi(xS,time,Global,Param)
          PsiCG[iD,iz] = PsiLoc
        end  
      end  
      @views Curl!(v1CG,v2CG,PsiCG,CG,Metric.dXdxI[:,:,:,:,:,iF],
        Metric.J[:,:,:,iF],Global.ThreadCache)
    else
      @. wCG = 0.0  
      @inbounds for iD = 1 : DoF
        @inbounds for iz=1:nz
          x1 = 0.5 * (X[iD,1,1,iz,iF] + X[iD,2,1,iz,iF])
          x2 = 0.5 * (X[iD,1,2,iz,iF] + X[iD,2,2,iz,iF])
          x3 = 0.5 * (X[iD,1,3,iz,iF] + X[iD,2,3,iz,iF])
          xS = SVector{3}(x1, x2 ,x3)
          uu = fVelu(xS,time,Phys,Global,Param)
          vv = fVelv(xS,time,Phys,Global,Param)
          v1CG[iD,iz] = uu
          v2CG[iD,iz] = vv
        end
        @inbounds for iz=1:nz-1
          x1 = 0.5 * (X[iD,2,1,iz,iF] + X[iD,1,1,iz+1,iF])
          x2 = 0.5 * (X[iD,2,2,iz,iF] + X[iD,1,2,iz+1,iF])
          x3 = 0.5 * (X[iD,2,3,iz,iF] + X[iD,1,3,iz+1,iF])
          xS = SVector{3}(x1, x2 ,x3)
          ww = fVelW(xS,time,Phys,Global,Param)
          wCG[iD,iz+1] = ww
        end  
      end  
    end
    @inbounds for iD = 1 : DoF
      ind = CG.Glob[iD,iF]
      @inbounds for iz=1:nz
        RhoCG[iD,iz] = U[iz,ind,RhoPos]
        @inbounds for iT = 1:NumTr
          TrCG[iD,iz,iT] = U[iz,ind,iT+NumV]
        end  
      end
    end
    @. FCG = 0.0

    @views DivRhoColumn!(FCG[:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Metric.dXdxI[:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))

#   Tracer transport
    @inbounds for iT = 1:NumTr
      @inbounds for iD = 1 : DoF
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          DivTrCG[iD,iz] = DivTr[iz,ind,iT] 
        end
      end
#     Hyperdiffusion, second Laplacian
      @views DivRhoGrad!(FCG[:,:,iT+NumV],DivTrCG,RhoCG,CG,
        Metric.dXdxI[:,:,:,:,:,iF],Metric.J[:,:,:,iF],Global.ThreadCache,
        Global.Model.HyperDDiv)
      if Global.Model.Upwind
        @views DivUpwindRhoTrColumn!(FCG[:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,iT],RhoCG,CG,
          Metric.dXdxI[:,:,:,:,:,iF],Metric.J[:,:,:,iF],
          Global.ThreadCache,Global.Model.HorLimit,Val(:VectorInvariant))
      else
        @views DivRhoTrColumn!(FCG[:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,iT],CG,
          Metric.dXdxI[:,:,:,:,:,],Global.ThreadCache,Val(:VectorInvariant))  
      end
    end


    
    if HorLimit
      @inbounds for iz=1:nz
        @inbounds for iT = 1:NumTr
          @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:],iT])
          @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:],iT])
          @views HorLimiter!(FCG[:,iz,iT+NumV],TrCG[:,iz,iT],RhoCG,RhoCG,dtau,
            Metric.J[:,:,iz,iF],CG.w,qMinS,qMaxS,Global.ThreadCache,iF)
        end
      end  
    end  

    @inbounds for iD = 1 : DoF
      ind = CG.Glob[iD,iF]
      @inbounds for iz=1:nz
        F[iz,ind,RhoPos] += FCG[iD,iz,RhoPos] / CG.M[iz,ind]
        @inbounds for iT = 1:NumTr
          F[iz,ind,iT+NumV] += FCG[iD,iz,iT+NumV] / CG.M[iz,ind]
        end
      end  
    end
  end  
  
  Parallels.ExchangeData3DSend(F,Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    if Param.StreamFun
      @inbounds for iD = 1 : DoF
        @inbounds for iz=1:nz
          x1 = 0.5 * (X[iD,1,1,iz,iF] + X[iD,2,1,iz,iF])
          x2 = 0.5 * (X[iD,1,2,iz,iF] + X[iD,2,2,iz,iF])
          x3 = 0.5 * (X[iD,1,3,iz,iF] + X[iD,2,3,iz,iF])
          xS = SVector{3}(x1, x2 ,x3)
          PsiLoc = fPsi(xS,time,Global,Param)
          PsiCG[iD,iz] = PsiLoc
        end 
      end  
      @views Curl!(v1CG,v2CG,PsiCG,CG,Metric.dXdxI[:,:,:,:,:,iF],
        Metric.J[:,:,:,iF],Global.ThreadCache)
    else
      @. wCG = 0.0  
      @inbounds for iD = 1 : DoF
        @inbounds for iz=1:nz
          x1 = 0.5 * (X[iD,1,1,iz,iF] + X[iD,2,1,iz,iF])
          x2 = 0.5 * (X[iD,1,2,iz,iF] + X[iD,2,2,iz,iF])
          x3 = 0.5 * (X[iD,1,3,iz,iF] + X[iD,2,3,iz,iF])
          xS = SVector{3}(x1, x2 ,x3)
          uu = Models.fVelu(xS,time,Phys,Global,Param)
          vv = Models.fVelv(xS,time,Phys,Global,Param)
          v1CG[iD,iz] = uu
          v2CG[iD,iz] = vv
        end 
        @inbounds for iz=1:nz-1
          x1 = 0.5 * (X[iD,2,1,iz,iF] + X[iD,1,1,iz+1,iF])
          x2 = 0.5 * (X[iD,2,2,iz,iF] + X[iD,1,2,iz+1,iF])
          x3 = 0.5 * (X[iD,2,3,iz,iF] + X[iD,1,3,iz+1,iF])
          xS = SVector{3}(x1, x2 ,x3)
          ww = Models.fVelW(xS,time,Phys,Global,Param)
          wCG[iD,iz+1] = ww
        end
      end  
    end 
    @inbounds for iD = 1 : DoF
      ind = CG.Glob[iD,iF]
      @inbounds for iz=1:nz
        RhoCG[iD,iz] = U[iz,ind,RhoPos]
        @inbounds for iT = 1:NumTr
          TrCG[iD,iz,iT] = U[iz,ind,iT+NumV]
        end  
      end
    end
    @. FCG = 0.0

    @views DivRhoColumn!(FCG[:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Metric.dXdxI[:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))

#   Tracer transport
    @inbounds for iT = 1 : NumTr
      @inbounds for iD = 1 : DoF
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          DivTrCG[iD,iz] = DivTr[iz,ind,iT] 
        end
      end
#     Hyperdiffusion, second Laplacian
      @views DivRhoGrad!(FCG[:,:,iT+NumV],DivTrCG,RhoCG,CG,
        Metric.dXdxI[:,:,:,:,:,iF],Metric.J[:,:,:,iF],Global.ThreadCache,
        Global.Model.HyperDDiv)
      if Global.Model.Upwind
        @views DivUpwindRhoTrColumn!(FCG[:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,iT],RhoCG,CG,
          Metric.dXdxI[:,:,:,:,:,iF],Metric.J[:,:,:,iF],
          Global.ThreadCache,Global.Model.HorLimit,Val(:VectorInvariant))
      else
        @views DivRhoTrColumn!(FCG[:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,iT],CG,
          Metric.dXdxI[:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))  
      end
    end

    if HorLimit
      @inbounds for iz=1:nz
        @inbounds for iT = 1:NumTr
          @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:],iT])
          @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:],iT])
          @views HorLimiter!(FCG[:,iz,iT+NumV],TrCG[:,iz,iT],RhoCG,RhoCG,dtau,
            Metric.J[:,:,iz,iF],CG.w,qMinS,qMaxS,Global.ThreadCache,iF)
        end
      end
    end  
      
    @inbounds for iD = 1 : DoF
      ind = CG.Glob[iD,iF]
      @inbounds for iz=1:nz
        F[iz,ind,RhoPos] += FCG[iD,iz,RhoPos] / CG.M[iz,ind]
        @inbounds for iT = 1:NumTr
          F[iz,ind,iT+NumV] += FCG[iD,iz,iT+NumV] / CG.M[iz,ind]
        end
      end
    end
  end  

  Parallels.ExchangeData3DRecv!(F,Exchange)

end

function FcnTracerConv!(F,U,time::Float64,CG,Metric,Global,Param)
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
  J = Metric.J
  X = Metric.X
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
      @views Curl!(v1CG,v2CG,PsiCG,CG,Metric.dXdxI[:,:,:,:,:,:,iF],
        Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @. wCG = 0.0
    else
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            x[1] = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x[2] = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x[3] = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            uu,vv = @gc_preserve Models.fVel(x,time,Global,Param)
            v1CG[iP,jP,iz] = uu
            v2CG[iP,jP,iz] = vv
          end
        end
      end
      @. wCG = 0.0
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz-1
            x[1] = 0.5 * (X[iP,jP,2,1,iz,iF] + X[iP,jP,1,1,iz+1,iF])
            x[2] = 0.5 * (X[iP,jP,2,2,iz,iF] + X[iP,jP,1,2,iz+1,iF])
            x[3] = 0.5 * (X[iP,jP,2,3,iz,iF] + X[iP,jP,1,3,iz+1,iF])
            ww = @gc_preserve Models.fVelw(x,time,Global,Param)
            wCG[iP,jP,iz+1] = uu
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
      Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))

    @inbounds for iT=1:NumTr
      @views DivConvRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],RhoCG,CG,
              Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.ThreadCache)
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
            x[1] = 0.5 * (X[iP,jP,1,1,iz,iF] + X[iP,jP,2,1,iz,iF])
            x[2] = 0.5 * (X[iP,jP,1,2,iz,iF] + X[iP,jP,2,2,iz,iF])
            x[3] = 0.5 * (X[iP,jP,1,3,iz,iF] + X[iP,jP,2,3,iz,iF])
            PsiLoc = @gc_preserve fPsi(x,time,Global,Param)
            PsiCG[iP,jP,iz] = PsiLoc
          end
        end
      end
      @views Curl!(v1CG,v2CG,PsiCG,CG,Metric.dXdxI[:,:,:,:,:,:,iF],
        Metric.J[:,:,:,:,iF],Global.ThreadCache)
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
      Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))

    @inbounds for iT = 1 : NumTr
      @views DivConvRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],RhoCG,CG,
        Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.ThreadCache)
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

function FcnTracerTest!(F,U,time,CG,Metric,Global,Param)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    NumV,
    NumTr) = Global.Model

  dtau = Global.TimeStepper.dtauStage
  HorLimit = Global.Model.HorLimit
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  J = Metric.J
  X = Metric.X
  Temp1 = Global.Cache.Temp1
  @. Temp1 = 0
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
  F .= 0.0


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
#     @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
#       dXdxI[:,:,:,:,:,:,iF],J[:,:,:,:,iF],Global.ThreadCache)
      @views DivRhoGradE!(DivCG,ThCG,RhoCG,CG.DS,CG.DW,
        Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF])
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
#     @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
#       Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @views DivRhoGradE!(DivCG,ThCG,RhoCG,CG.DS,CG.DW,
        Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF])
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

  @views ExchangeData3DRecv!(Temp1[:,:,1:NumV+NumTr],Global.Exchange)

  @views Velocity!(U[:,:,uPos],U[:,:,vPos],U[:,:,wPos],time,CG.Glob,X,Param)

  @inbounds for iF in Global.Grid.BoundaryFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          wCG[iP,jP,iz+1] = U[iz,ind,wPos]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,iT+NumV]
          end  
        end
      end
    end
    @. FCG = 0.0

#   @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
#     Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))
### @views DivRhoColumnE!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG.DS,
###   Metric.dXdxI[:,:,:,:,:,:,iF])

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
#     @views DivRhoGrad!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG,
#       Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.ThreadCache,
#       Global.Model.HyperDDiv)
      @views DivRhoGradE!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG.DS,CG.DW,
        Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.Model.HyperDDiv)
      if Global.Model.Upwind
        @views DivUpwindRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],RhoCG,CG,
          Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],
          Global.ThreadCache,Global.Model.HorLimit,Val(:VectorInvariant))
      else
###     @views DivRhoTrColumnE!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG.DS,
###     Metric.dXdxI[:,:,:,:,:,:,iF])
#       @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
#         Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))  
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
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          wCG[iP,jP,iz+1] = U[iz,ind,wPos]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,iT+NumV]
          end  
        end
      end
    end
    @. FCG = 0.0

#   @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
#     Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))
### @views DivRhoColumnE!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG.DS,
###   Metric.dXdxI[:,:,:,:,:,:,iF])

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
      @views DivRhoGrad!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG,
        Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.ThreadCache,
        Global.Model.HyperDDiv)
###   @views DivRhoGradE!(FCG[:,:,:,iT+NumV],DivTrCG,RhoCG,CG.DS,CG.DW,
###     Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],Global.Model.HyperDDiv)
      if Global.Model.Upwind
        @views DivUpwindRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],RhoCG,CG,
          Metric.dXdxI[:,:,:,:,:,:,iF],Metric.J[:,:,:,:,iF],
          Global.ThreadCache,Global.Model.HorLimit,Val(:VectorInvariant))
      else
#       @views DivRhoTrColumnE!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG.DS,
#         Metric.dXdxI[:,:,:,:,:,:,iF])
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
          Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache,Val(:VectorInvariant))  
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

function Velocity!(u,v,w,time,Glob,X,Param)
  Profile = RotationalCart(RotationalCartExample())
  NumFaces = size(X,6)
  Nz = size(X,5)
  N = size(X,1)
  for iF = 1 : NumFaces
    for iz = 1 : Nz
      for j = 1 : N
        for i = 1 : N 
          ind = Glob[i,j,iF]
          x1 = 0.5 * (X[i,j,1,1,iz,iF] + X[i,j,2,1,iz,iF])
          x2 = 0.5 * (X[i,j,1,2,iz,iF] + X[i,j,2,2,iz,iF])
          x3 = 0.5 * (X[i,j,1,3,iz,iF] + X[i,j,2,3,iz,iF])
          xS = SVector{3}(x1, x2 ,x3)
          _,u[iz,ind],v[iz,ind],_ = Profile(xS,time,Param)
        end
      end
    end
    for iz = 1 : Nz - 1
      for j = 1 : N
        for i = 1 : N 
          ind = Glob[i,j,iF]
          x1 = 0.5 * (X[i,j,2,1,iz,iF] + X[i,j,1,1,iz+1,iF])
          x2 = 0.5 * (X[i,j,2,2,iz,iF] + X[i,j,1,2,iz+1,iF])
          x3 = 0.5 * (X[i,j,2,3,iz,iF] + X[i,j,1,3,iz+1,iF])
          xS = SVector{3}(x1, x2 ,x3)
          _,_,_,w[iz,ind] = Profile(xS,time,Param)
        end
      end
    end
  end
end
