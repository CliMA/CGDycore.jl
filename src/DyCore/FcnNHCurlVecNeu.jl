function FcnNHCurlVecINeu!(F,U,CG,Global,Param)

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
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  @views JJ = Global.Cache.Temp1[:,:,5+NumTr+1]
  @views JRho = Global.Cache.Temp1[:,:,5+NumTr+2]
  @views JRhoF = Global.Cache.Temp1[:,:,5+NumTr+3]
  FCG=Global.Cache.FCC
  FwCG=Global.Cache.FwCC
  Rot1CG=Global.Cache.Rot1C
  Rot2CG=Global.Cache.Rot2C
  Grad1CG=Global.Cache.Grad1C
  Grad2CG=Global.Cache.Grad2C
  DivCG=Global.Cache.DivC
  DivTrCG=Global.Cache.DivC
  @views RhoCG = Global.Cache.RhoCG[:,:,:]
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  zPG = Global.Cache.zPG
  pBGrdCG = Global.Cache.pBGrdCG
  RhoBGrdCG = Global.Cache.RhoBGrdCG
  @views ThCG = Global.Cache.ThCG[:,:,:]
  @views TrCG = Global.Cache.TrCG[:,:,:,:]
  KE = Global.Cache.KE
  Pres = Global.Cache.Pres
  PresG = Global.Cache.PresG
  Temp = Global.Cache.Temp
  uStar = Global.Cache.uStar
  JC = Global.Metric.JC
  KV = Global.Cache.KV
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  DivTr .= 0.0
  F .= 0.0
  PresG .= 0.0
  JJ .= 0.0
  JRho .= 0.0
  JRhoF .= 0.0
  # Hyperdiffusion 
  @inbounds for iF in Global.Grid.BoundaryFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
        end
      end
    end
      

    @views RotCurl!(Rot1CG,Rot2CG,v1CG,v2CG,CG,
     Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    @views GradDiv!(Rot1CG,Rot2CG,v1CG,v2CG,CG,
     Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
     Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] 
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] 
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] 
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] 
          Div[iz,ind] += DivCG[iP,jP,iz] 
          JJ[iz,ind] += J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]
          JRho[iz,ind] += (J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]) * RhoCG[iP,jP,iz]
        end
        @inbounds for iz=1:nz-1
          JRhoF[iz,ind] += (J[iP,jP,2,iz,iF] * RhoCG[iP,jP,iz] + J[iP,jP,1,iz+1,iF] * RhoCG[iP,jP,iz+1])
        end
      end
    end
    @inbounds for iT=1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            ThCG[iP,jP,iz] = U[iz,ind,NumV+iT]
          end
        end
      end
      @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTr[iz,ind,iT] += DivCG[iP,jP,iz]
          end
        end
      end
    end
  end

  ExchangeData3DSend(Temp1,PresG,Global.Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
        end
      end
    end
      

    @views RotCurl!(Rot1CG,Rot2CG,v1CG,v2CG,CG,
     Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    @views GradDiv!(Rot1CG,Rot2CG,v1CG,v2CG,CG,
     Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
     Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] 
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] 
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] 
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] 
          Div[iz,ind] += DivCG[iP,jP,iz] 
          JJ[iz,ind] += J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]
          JRho[iz,ind] += (J[iP,jP,1,iz,iF] + J[iP,jP,2,iz,iF]) * RhoCG[iP,jP,iz]
        end
        @inbounds for iz=1:nz-1
          JRhoF[iz,ind] += (J[iP,jP,2,iz,iF] * RhoCG[iP,jP,iz] + J[iP,jP,1,iz+1,iF] * RhoCG[iP,jP,iz+1])
        end
      end
    end
    @inbounds for iT=1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            ThCG[iP,jP,iz] = U[iz,ind,NumV+iT]
          end
        end
      end
      @views DivRhoGrad!(DivCG,ThCG,RhoCG,CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivTr[iz,ind,iT] += DivCG[iP,jP,iz] 
          end
        end
      end
    end
  end

  ExchangeData3DRecv!(Temp1,PresG,Global.Exchange)

  @inbounds for iF in Global.Grid.BoundaryFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          wCG[iP,jP,iz+1] = U[iz,ind,wPos]
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
          Rot1CG[iP,jP,iz] = Rot1[iz,ind] / JJ[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind] / JJ[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind] / JJ[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind] / JJ[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind] / JJ[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          pBGrdCG[iP,jP,iz] = Global.pBGrd[iz,ind]
          RhoBGrdCG[iP,jP,iz] = Global.RhoBGrd[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT] / JJ[iz,ind]
          end  
        end
      end
    end
    @. FCG = 0.0
    @. FwCG = 0.0

#   Hyperdiffusion Part 2
    @views RotCurl!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],Rot1CG,Rot2CG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
      Global.Model.HyperDCurl)
    @views GradDiv!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],Grad1CG,Grad2CG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
      Global.Model.HyperDGrad)
    @views DivRhoGrad!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
      Global.Model.HyperDDiv)

#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    @views BoundaryW!(wCG[:,:,1],v1CG[:,:,1],v2CG[:,:,1],CG,
      Global.Metric.dXdxI[:,:,1,1,:,:,iF])
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zPG,Global)
#   Temperature

    if Global.Model.VerticalDiffusion || Global.Model.SurfaceFlux
#     uStar
      @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,iF)

#     Vertical Diffusion coefficient    
      if Global.Model.VerticalDiffusion
        eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,Param,iF)
      end   
    end   

#   @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);
    @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

    if Global.Model.RefProfile
  #   Global.Pres=Pressure(ThCG,RhoCG,KE,Global)-Global.pBGrd;
  #   FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Global.Pres,CG,Global);
  #   FCG[:,:,:,uPos]=FCG[:,:,:,uPos]./RhoCG;
  #   FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
  #   if Global.Buoyancy
  #     RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
  #     FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
  #       Global.Grav*Global.JF[:,:,:,2:nz].*(RhoF-Global.RhoBGrdF)) ./ RhoF;
  #   end
    else
      @views @. pBGrdCG = Pres[:,:,:,iF] - pBGrdCG  
      @views GradColumn!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],Pres[:,:,:,iF],CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
      if Global.Model.Buoyancy
        @views Buoyancy!(FwCG,RhoCG,Global.Metric.J[:,:,:,:,iF],Phys)  
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views RhoGradColumn!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],KE,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    @views CurlColumn!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],
      v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "TotalEnergy"
      @views @. ThCG = ThCG + Pres[:,:,:,iF]  
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    elseif Global.Model.Thermo == "InternalEnergy" 
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
      @views SourceIntEnergy!(FCG[:,:,:,ThPos],Pres[:,:,:,iF],v1CG,v2CG,wCG,CG,Global,iF)
    else
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
#       @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
        @views DivRhoTrColumn!(FCG[:,:,:,ThPos],v1CG,v2CG,wCG,ThCG,CG,
          Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)
      end
      if Global.Model.VerticalDiffusion
        @views VerticalDiffusionScalar!(FCG[:,:,:,ThPos],ThCG,RhoCG,KV,CG,Global,iF)
      end  
    end  
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
#       @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
          Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)
      end
      if Global.Model.VerticalDiffusion
        @views VerticalDiffusionScalar!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],RhoCG,KV,CG,Global,iF)
      end  
    end
    if Global.Model.SurfaceFlux
      @views BoundaryFluxScalar!(FCG[:,:,1,:],ThCG[:,:,1],RhoCG[:,:,1],TrCG[:,:,1,:],CG,Global,Param,iF)
    end  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] 
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] 
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] 
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] 
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FwCG[iP,jP,iz+1] 
        end
      end  
    end
  end  

  ExchangeData3DSend(F,PresG,Global.Exchange)

  @inbounds for iF in Global.Grid.InteriorFaces
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          wCG[iP,jP,iz+1] = U[iz,ind,wPos]
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
          Rot1CG[iP,jP,iz] = Rot1[iz,ind] / JJ[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind] / JJ[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind] / JJ[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind] / JJ[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind] / JJ[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0
    @. FwCG = 0.0

    #   Hyperdiffusion Part 2
    @views RotCurl!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],Rot1CG,Rot2CG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
      Global.Model.HyperDCurl)
    @views GradDiv!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],Grad1CG,Grad2CG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
      Global.Model.HyperDGrad)
    @views DivRhoGrad!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache,
      Global.Model.HyperDDiv)

#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    @views BoundaryW!(wCG[:,:,1],v1CG[:,:,1],v2CG[:,:,1],CG,
      Global.Metric.dXdxI[:,:,1,1,:,:,iF])
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zPG,Global)
#   Temperature

    if Global.Model.VerticalDiffusion || Global.Model.SurfaceFlux
#     uStar
      @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,iF)

#     Vertical Diffusion coefficient    
      if Global.Model.VerticalDiffusion
        eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,Param,iF)
      end   
    end   

    @views DivRhoColumn!(FCG[:,:,:,RhoPos],v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

    if Global.Model.RefProfile
  #   Global.Pres=Pressure(ThCG,RhoCG,KE,Global)-Global.pBGrd;
  #   FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Global.Pres,CG,Global);
  #   FCG[:,:,:,uPos]=FCG[:,:,:,uPos]./RhoCG;
  #   FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
  #   if Global.Buoyancy
  #     RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
  #     FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
  #       Global.Grav*Global.JF[:,:,:,2:nz].*(RhoF-Global.RhoBGrdF)) ./ RhoF;
  #   end
    else
      @views GradColumn!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],Pres[:,:,:,iF],CG,
        Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)

      if Global.Model.Buoyancy
        @views Buoyancy!(FwCG,RhoCG,Global.Metric.J[:,:,:,:,iF],Phys)  
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views RhoGradColumn!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],KE,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)
    @views CurlColumn!(FCG[:,:,:,uPos],FCG[:,:,:,vPos],FwCG[:,:,:],
      v1CG,v2CG,wCG,RhoCG,CG,
      Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "TotalEnergy"
      @views @. ThCG = ThCG + Pres[:,:,:,iF]  
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    elseif Global.Model.Thermo == "InternalEnergy" 
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
      @views SourceIntEnergy!(FCG[:,:,:,ThPos],Pres[:,:,:,iF],v1CG,v2CG,wCG,CG,Global,iF)
    else
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views DivRhoTrColumn!(FCG[:,:,:,ThPos],v1CG,v2CG,wCG,ThCG,CG,
          Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)

      end
      if Global.Model.VerticalDiffusion
        @views VerticalDiffusionScalar!(FCG[:,:,:,ThPos],ThCG,RhoCG,KV,CG,Global,iF)
      end  
    end  
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
#       @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
        @views DivRhoTrColumn!(FCG[:,:,:,iT+NumV],v1CG,v2CG,wCG,TrCG[:,:,:,iT],CG,
          Global.Metric.dXdxI[:,:,:,:,:,:,iF],Global.ThreadCache)
      end
      if Global.Model.VerticalDiffusion
        @views VerticalDiffusionScalar!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],RhoCG,KV,CG,Global,iF)
      end  
    end
    if Global.Model.SurfaceFlux
      @views BoundaryFluxScalar!(FCG[:,:,1,:],ThCG[:,:,1],RhoCG[:,:,1],TrCG[:,:,1,:],CG,Global,Param,iF)
    end  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] 
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] 
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] 
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FwCG[iP,jP,iz+1] 
        end
      end  
    end
  end  
  ExchangeData3DRecv!(F,PresG,Global.Exchange)
  @views @. F[:,:,RhoPos] /= JJ
  @views @. F[:,:,ThPos] /= JJ
  @inbounds for iT = 1:NumTr
    @views @.  F[:,:,iT+NumV] /= JJ
  end
  @views @. F[:,:,uPos] /= JRho
  @views @. F[:,:,vPos] /= JRho
  @views @. F[1:nz-1,:,wPos] /= JRhoF[1:nz-1,:]

  if Global.Model.Damping
    @inbounds for iG=1:CG.NumG
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
  end
  if Global.Model.Source
    @inbounds for iG=1:CG.NumG
      @views Source!(F[:,iG,:],U[:,iG,:],PresG[:,iG],CG,Global,Param,iG)
    end
  end
  if Global.Model.Microphysics
    @inbounds for iG=1:CG.NumG
      @views SourceMicroPhysics(F[:,iG,:],U[:,iG,:],PresG[:,iG],CG,Global,iG)
    end  
  end
end

