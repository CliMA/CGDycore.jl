  function FcnTracerColor!(F,Tr,U,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TTrCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumTr,
    NumV) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  DivTr = Global.Cache.DivTr
  DivTr .= 0.0
  F .= 0.0
  # Hyperdiffusion 
  @inbounds for color in Global.Grid.colors
    @inbounds for iF in color
      TrCG = TTrCG[Threads.threadid()]
      DivTrCG = TDivTrCG[Threads.threadid()]
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind=CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            @inbounds for iT=1:NumTr
              TrCG[iz,ind,iT] = Tr[iz,ind,iT]   
            end  
          end
        end
      end
      FDivGrad2VecDSS!(DivTrCG,TrCG,RhoCG,CG,Global,iF)  
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            @inbounds for iT=1:NumTr
              DivTr[iz,ind,iT] += DivCG[iP,jP,iz,iT] / CG.M[iz,ind]
            end  
          end
        end
      end
    end
  end

  @inbounds for color in Global.Grid.colors
#   Threads.@threads for iF in color
    @inbounds for iF in color
      RhoCG = TRhoCG[Threads.threadid()]
      v1CG = Tv1CG[Threads.threadid()]
      v2CG = Tv2CG[Threads.threadid()]
      wCG  = TwCG[Threads.threadid()]
      wCCG  = TwCCG[Threads.threadid()]
      DivTrCG = TDivTrCG[Threads.threadid()]
      FCG = TFCG[Threads.threadid()]

      @. FCG = 0.0  
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
            v1CG[iP,jP,iz] = U[iz,ind,uPos]
            v2CG[iP,jP,iz] = U[iz,ind,vPos]
            wCG[iP,jP,iz+1] = U[iz,ind,wPos]
            @inbounds for iT=1:NumTr
              DivTrCG[iz,ind,iz,iT] = DivTR[iz,ind,iT] / CG.M[iz,ind]
            end  
            
          end
        end
      end
      @views FDivRhoGrad2Vec!(FCG,DivTrCG,RhoCG,CG,Global,iF)
      @views BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
      @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])

#     Divergence of Tracer Variable
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG,TrCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG,TrCG,v1CG,v2CG,wCG,CG,Global,iF);
      end  
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            @inbouns iT=1:NumTr
            F[iz,ind,iT] += FCG[iP,jP,iz,iT] / CG.M[iz,ind]
          end
        end  
      end
    end  
  end
end

