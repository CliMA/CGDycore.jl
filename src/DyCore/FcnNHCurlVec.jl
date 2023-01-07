function FcnNHCurlVec!(F,U,CG,Global,Param)

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  Temp1 = Global.Cache.Temp1
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  FCG=Global.Cache.FCC
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
  @views ThCG = Global.Cache.ThCG[:,:,:]
  @views TrCG = Global.Cache.TrCG[:,:,:,:]
  KE = Global.Cache.KE
  Pres = Global.Cache.Pres
  Temp = Global.Cache.Temp
  uStar = Global.Cache.uStar
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  DivTr .= 0.0
  F .= 0.0
  KV = Global.Cache.KV
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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

  ExchangeData3D!(Temp1,Global.Exchange)

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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zP[:,:,:,iF],Global)
#   Temperature

    if Global.Model.VerticalDiffusion || Global.Model.SurfaceFlux
#     uStar
      @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,Param,iF)

#     Vertical Diffusion coefficient    
      if Global.Model.VerticalDiffusion
        eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,Param,iF)
      end   
    end   

    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "Energy"
    else
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
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
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  

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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zP[:,:,:,iF],Global)
#   Temperature

    if Global.Model.VerticalDiffusion || Global.Model.SurfaceFlux
#     uStar
      @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,iF)

#     Vertical Diffusion coefficient    
      if Global.Model.VerticalDiffusion
        eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,Param,iF)
      end   
    end   

    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "Energy"
    else
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
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
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  

  ExchangeData3D!(F,Global.Exchange)

  if Global.Model.Damping
    @inbounds for iG=1:CG.NumG
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
  end
  if Global.Model.Source
    @inbounds for iG=1:CG.NumG
      @views Source!(F[:,iG,:],U[:,iG,:],CG,Global,Param,iG)
    end
  end
  if Global.Model.Microphysics
    @inbounds for iG=1:CG.NumG
      @views SourceMicroPhysics(F[:,iG,:],U[:,iG,:],CG,Global,iG)
    end  
  end
end

function FcnNHCurlVecI!(F,U,CG,Global,Param)

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  FCG=Global.Cache.FCC
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          pBGrdCG[iP,jP,iz] = Global.pBGrd[iz,ind]
          RhoBGrdCG[iP,jP,iz] = Global.RhoBGrd[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
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

    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      #@views FGrad3RhoVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      @views FGrad3RhoVec!(FCG,pBGrdCG,RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              Buo = (RhoCG[i,j,iz] + RhoCG[i,j,iz+1] - RhoBGrdCG[i,j,iz] - RhoBGrdCG[i,j,iz+1]) /
                (RhoCG[i,j,iz] + RhoCG[i,j,iz+1]) 
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF] * Buo
            end
          end  
        end
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

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
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
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
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
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

    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

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
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
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
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  
  ExchangeData3DRecv!(F,PresG,Global.Exchange)

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

function FcnNHCurlExp3DVecI!(F,U,CG,Global,Param)

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  FCG=Global.Cache.FCC
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure

    if Global.Model.VerticalDiffusion || Global.Model.SurfaceFlux
#     uStar
      @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,iF)

#     Vertical Diffusion coefficient    
      if Global.Model.VerticalDiffusion
        eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,Param,iF)
      end   
    end   

#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "TotalEnergy"
    elseif Global.Model.Thermo == "InternalEnergy" 
    else
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
    end
    if Global.Model.SurfaceFlux
      @views BoundaryFluxScalar!(FCG[:,:,1,:],ThCG[:,:,1],RhoCG[:,:,1],TrCG[:,:,1,:],CG,Global,Param,iF)
    end  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))

    if Global.Model.VerticalDiffusion || Global.Model.SurfaceFlux
#     uStar
      @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,iF)

#     Vertical Diffusion coefficient    
      if Global.Model.VerticalDiffusion
        eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,Param,iF)
      end   
    end   

#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "TotalEnergy"
    elseif Global.Model.Thermo == "InternalEnergy" 
    else
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
    end
    if Global.Model.SurfaceFlux
      @views BoundaryFluxScalar!(FCG[:,:,1,:],ThCG[:,:,1],RhoCG[:,:,1],TrCG[:,:,1,:],CG,Global,Param,iF)
    end  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  
  ExchangeData3DRecv!(F,PresG,Global.Exchange)

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

function FcnNHCurlImp3DVecI!(F,U,CG,Global,Param)

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  FCG=Global.Cache.FCC
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
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zPG,Global)
#   Temperature

    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end

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
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    end  
#   Tracer transport
    @inbounds for iT = 1:NumTr
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
      end
      if Global.Model.VerticalDiffusion
        @views VerticalDiffusionScalar!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],RhoCG,KV,CG,Global,iF)
      end  
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
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
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zPG,Global)
#   Temperature

    @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end

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
        @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    end  
#   Tracer transport
    @inbounds for iT = 1:NumTr
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
      end
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end
  ExchangeData3DRecv!(F,PresG,Global.Exchange)
end

function FcnNHCurlExp1DVecI!(F,U,CG,Global,Param)

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  FCG=Global.Cache.FCC
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
      

    @views FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

    @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
      @views FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)
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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
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

    @views FDiv3ExpVec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoExpVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

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
      @views FDiv3ExpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
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
      @views FDiv3ExpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
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
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Hyperdiffusion Part 2    
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)


#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
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

    @views FDiv3ExpVec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoExpVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG,KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

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
      @views FDiv3ExpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
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
      @views FDiv3ExpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
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
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  
  ExchangeData3DRecv!(F,PresG,Global.Exchange)

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

function FcnNHCurlImp1DVecI!(F,U,CG,Global,Param)

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV,
    NumTr) = Global.Model
  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  zP = Global.Metric.zP
  Temp1 = Global.Cache.Temp1
  @views Rot1 = Global.Cache.Temp1[:,:,1]
  @views Rot2 = Global.Cache.Temp1[:,:,2]
  @views Grad1 = Global.Cache.Temp1[:,:,3]
  @views Grad2 = Global.Cache.Temp1[:,:,4]
  @views Div = Global.Cache.Temp1[:,:,5]
  @views DivTr = Global.Cache.Temp1[:,:,5+1:5+NumTr]
  FCG=Global.Cache.FCC
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
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zPG,Global)
#   Temperature

    @views FDiv3ImpVec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoImpVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end

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
        @views FDiv3UpwindImpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3ImpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    end  
#   Tracer transport
    @inbounds for iT = 1:NumTr
      if Global.Model.Upwind
        @views FDiv3UpwindImpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3ImpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
      end
      if Global.Model.VerticalDiffusion
        @views VerticalDiffusionScalar!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],RhoCG,KV,CG,Global,iF)
      end  
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
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
          zPG[iP,jP,iz] = zP[iz,ind]
          @inbounds for iT = 1:NumTr
            TrCG[iP,jP,iz,iT] = U[iz,ind,NumV+iT]
          end  
        end
      end
    end
    @. FCG = 0.0

#   Diagnostic values
#   Boundary value for vertical velocity and cell center   
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#   Pressure
    @views Pressure!(Pres[:,:,:,iF],ThCG,RhoCG,TrCG,KE,zPG,Global)
#   Temperature

    @views FDiv3ImpVec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

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
      @views FGrad3RhoImpVec!(FCG,Pres[:,:,:,iF],RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end

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
        @views FDiv3UpwindImpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3ImpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    end  
#   Tracer transport
    @inbounds for iT = 1:NumTr
      if Global.Model.Upwind
        @views FDiv3UpwindImpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3ImpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
      end
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          PresG[iz,ind,RhoPos] += Pres[iP,jP,iz,iF] * JC[iP,jP,iz,iF] / CG.M[iz,ind]
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          @inbounds for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end
  ExchangeData3DRecv!(F,PresG,Global.Exchange)
end

function FcnNHCurlImp1DGlobalVecI!(F,U,CG,Global,Param)

  (;  RhoPos,
      uPos,
      vPos,
      wPos,
      ThPos,
      NumV,
      NumTr) = Global.Model
  PresG = Global.Cache.PresG
  zP = Global.Metric.zP
  Grav=Global.Phys.Grav    
  nz = Global.Grid.nz
  @. F = 0.0
  @inbounds for iG=1:CG.NumG
    @views Pressure!(PresG[:,iG],U[:,iG,ThPos],U[:,iG,RhoPos],U[:,iG,NumV+1:NumV+NumTr],
      U[:,iG,ThPos],zP[:,iG],Global)
    @views FDiv3ImpGlobalVec!(F[:,iG,RhoPos],U[:,iG,RhoPos],U[:,iG,wPos],Global,iG)
    if Global.Model.Upwind
      @views FDiv3UpwindImpGlobalVec!(F[:,iG,ThPos],U[:,iG,ThPos],U[:,iG,wPos],U[:,iG,RhoPos],Global,iG)
    else
      @views FDiv3ImpGlobalVec!(F[:,iG,ThPos],U[:,iG,ThPos],U[:,iG,wPos],Global,iG)
    end
    @views FGrad3RhoImpGlobalVec!(F[:,iG,wPos],PresG[:,iG],U[:,iG,RhoPos],Global,iG)
    @views @. F[1:nz-1,iG,wPos] = F[1:nz-1,iG,wPos] - Grav  
#   Tracer transport
    @inbounds for iT = 1:NumTr
      if Global.Model.Upwind
        @views FDiv3UpwindImpGlobalVec!(F[:,iG,iT+NumV],U[:,iG,iT+NumV],U[:,iG,wPos],U[:,iG,RhoPos],Global,iG)
      else
        @views FDiv3ImpGlobalVec!(F[:,iG,iT+NumV],U[:,iG,iT+NumV],U[:,iG,wPos],Global,iG)
      end
    end
  end  

#   Divergence of Thermodynamic Variable
#   if Global.Model.Thermo == "TotalEnergy"
#     @views @. ThCG = ThCG + Pres[:,:,:,iF]  
#     if Global.Model.Upwind
#       @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
#     else
#       @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
#     end
#   elseif Global.Model.Thermo == "InternalEnergy" 
#     if Global.Model.Upwind
#       @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
#     else
#       @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
#     end
#     @views SourceIntEnergy!(FCG[:,:,:,ThPos],Pres[:,:,:,iF],v1CG,v2CG,wCG,CG,Global,iF)
#   else
#     if Global.Model.Upwind
#       @views FDiv3UpwindImpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
#     else
#       @views FDiv3ImpVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
#     end
#   end  
#   Tracer transport
#   @inbounds for iT = 1:NumTr
#     if Global.Model.Upwind
#       @views FDiv3UpwindImpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
#     else
#       @views FDiv3ImpVec!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],v1CG,v2CG,wCG,CG,Global,iF);
#     end
#   end

end

function FcnNHCurlVecIZero!(F,U,CG,Global,Param)
  @. F = 0.0
end  
