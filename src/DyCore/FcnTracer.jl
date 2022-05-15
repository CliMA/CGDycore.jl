function FcnTracer!(F,U,time,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    NumV,
    NumTr) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  DivTr = Global.Cache.DivTr
  DivCG=Global.Cache.DivC
  FCG=Global.Cache.FCC
  DivTrCG=Global.Cache.DivC
  RhoCG = Global.Cache.RhoCG
  ThCG = Global.Cache.ThCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  TrCG = Global.Cache.TrCG
  DivTr .= 0.0
  F .= 0.0
  if Global.Model.Param.TimeDependent
    @views ProjectVec!(U[:,:,uPos],U[:,:,vPos],fVel,time,CG,Global)
    @views ProjectW!(U[:,:,wPos],fVelW,time,CG,Global)
    @views @. U[:,CG.Boundary,NumV+1:end] = 0.0
    @views @. U[:,CG.Boundary,RhoPos] = 1.0
  end
  # Hyperdiffusion 
  @inbounds for iF = 1:NF
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

  @inbounds for iF = 1:NF
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
        @views FDiv3Vec!(FCG[:,:,:,iT],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
end

function ComputeVelocity!(U,CG,Global,time)
  Model=Global.Model
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Global)
end
