  function FcnNHCurlVecColor!(F,U,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  Rot1 = Global.Cache.Rot1
  Rot2 = Global.Cache.Rot2
  Grad1 = Global.Cache.Grad1
  Grad2 = Global.Cache.Grad2
  Div = Global.Cache.Div
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  F .= 0.0
  # Hyperdiffusion 
  @inbounds for color in Global.Grid.colors
    @inbounds for iF in color
      RhoCG = TRhoCG[Threads.threadid()]
      v1CG = Tv1CG[Threads.threadid()]
      v2CG = Tv2CG[Threads.threadid()]
      ThCG = TThCG[Threads.threadid()]
      Rot1CG = TCacheCC1[Threads.threadid()]
      Rot2CG = TCacheCC2[Threads.threadid()]
      Grad1CG = TCacheCC3[Threads.threadid()]
      Grad2CG = TCacheCC4[Threads.threadid()]
      DivCG = TCacheCC5[Threads.threadid()]
      iG=0
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          @inbounds for iz=1:nz
            RhoCG[iP,jP,iz] = U[iz,CG.Glob[iP,jP,iF],RhoPos]
            v1CG[iP,jP,iz] = U[iz,CG.Glob[iP,jP,iF],uPos]
            v2CG[iP,jP,iz] = U[iz,CG.Glob[iP,jP,iF],vPos]
            ThCG[iP,jP,iz] = U[iz,CG.Glob[iP,jP,iF],ThPos]
          end
        end
      end

      FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
      FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

#     FTempW=FDivGrad2Vec(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,Global);
      FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  
      iG=0
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
      ThCG = TThCG[Threads.threadid()]
      Rot1CG = TCacheCC1[Threads.threadid()]
      Rot2CG = TCacheCC2[Threads.threadid()]
      Grad1CG = TCacheCC3[Threads.threadid()]
      Grad2CG = TCacheCC4[Threads.threadid()]
      DivCG = TCacheCC5[Threads.threadid()]
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
            ThCG[iP,jP,iz] = U[iz,ind,ThPos]
            Rot1CG[iP,jP,iz] = Rot1[iz,ind]
            Rot2CG[iP,jP,iz] = Rot2[iz,ind]
            Grad1CG[iP,jP,iz] = Grad1[iz,ind]
            Grad2CG[iP,jP,iz] = Grad2[iz,ind]
            DivCG[iP,jP,iz] = Div[iz,ind]
          end
        end
      end
      @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
      @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
      @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)
      BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
      @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])

      KE = TCacheCC1[Threads.threadid()]
      @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + wCCG*wCCG);

      @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

      if Global.Model.RefProfile
  #     Global.Pres=Pressure(ThCG,RhoCG,KE,Global)-Global.pBGrd;
  #     FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Global.Pres,CG,Global);
  #     FCG[:,:,:,:,uPos]=FCG[:,:,:,:,uPos]./RhoCG;
  #     FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
  #     if Global.Buoyancy
  #       RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
  #       FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
  #         Global.Grav*Global.JF[:,:,:,2:nz].*(RhoF-Global.RhoBGrdF)) ./ RhoF;
  # end
      else
        Pres = TCacheCC2[Threads.threadid()]
        Pressure!(Pres,ThCG,RhoCG,KE,Global);
        FGrad3RhoVec!(FCG,Pres,RhoCG,CG,Global,iF)
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
#     3-dim Curl and Grad of kinetic Energy
      FGrad3Vec!(FCG,KE,CG,Global,iF)
      FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

#     Divergence of Thermodynamic Variable
      if Global.Model.Thermo == "Energy"
      else
        if Global.Model.Upwind
          @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
        else
          @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
        end
      end  
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
            F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
            F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
            F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          end
          @inbounds for iz=1:nz-1
            F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
          end
        end  
      end
    end  
  end

  @inbounds for iG=1:CG.NumG
    if Global.Model.Damping
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
    if Global.Model.Source
      @views Source!(F[:,iG,:],U[:,iG,:],CG,Global);
    end
  end
end

function FcnNHCurlVecThreads!(F,U,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  Rot1 = Global.Cache.Rot1
  Rot2 = Global.Cache.Rot2
  Grad1 = Global.Cache.Grad1
  Grad2 = Global.Cache.Grad2
  Div = Global.Cache.Div
  FCG=Global.Cache.FCG
  Rot1CG=Global.Cache.Rot1CG
  Rot2CG=Global.Cache.Rot2CG
  Grad1CG=Global.Cache.Grad1CG
  Grad2CG=Global.Cache.Grad2CG
  DivCG=Global.Cache.DivCG
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  F .= 0.0
  FCG .= 0.0
  # Hyperdiffusion 
# Threads.@threads for iF = 1:NF
  @inbounds for iF = 1:NF
    RhoCG = TRhoCG[Threads.threadid()]
    v1CG = Tv1CG[Threads.threadid()]
    v2CG = Tv2CG[Threads.threadid()]
    ThCG = TThCG[Threads.threadid()]

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

    @views FRotCurl2VecDSS!(Rot1CG[:,:,:,iF],Rot2CG[:,:,:,iF],v1CG,v2CG,CG,Global,iF)
    @views FGradDiv2VecDSS!(Grad1CG[:,:,:,iF],Grad2CG[:,:,:,iF],v1CG,v2CG,CG,Global,iF)

    #FTempW=FDivGrad2Vec(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,Global);
    @views FDivGrad2VecDSS!(DivCG[:,:,:,iF],ThCG,RhoCG,CG,Global,iF)  
  end    

  @inbounds for iF = 1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          Rot1[iz,ind] += Rot1CG[iP,jP,iz,iF] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz,iF] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz,iF] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz,iF] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz,iF] / CG.M[iz,ind]
        end
      end
    end
  end

  #Threads.@threads for iF = 1:NF
  for iF = 1:NF
    RhoCG = TRhoCG[Threads.threadid()]
    v1CG = Tv1CG[Threads.threadid()]
    v2CG = Tv2CG[Threads.threadid()]
    wCG  = TwCG[Threads.threadid()]
    wCCG  = TwCCG[Threads.threadid()]
    ThCG = TThCG[Threads.threadid()]
    Rot1C = TCacheCC1[Threads.threadid()]
    Rot2C = TCacheCC2[Threads.threadid()]
    Grad1C = TCacheCC3[Threads.threadid()]
    Grad2C = TCacheCC4[Threads.threadid()]
    DivC = TCacheCC5[Threads.threadid()]

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
          v1CG[iP,jP,iz] = U[iz,ind,uPos]
          v2CG[iP,jP,iz] = U[iz,ind,vPos]
          wCG[iP,jP,iz+1] = U[iz,ind,wPos]
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
          Rot1C[iP,jP,iz] = Rot1[iz,ind]
          Rot2C[iP,jP,iz] = Rot2[iz,ind]
          Grad1C[iP,jP,iz] = Grad1[iz,ind]
          Grad2C[iP,jP,iz] = Grad2[iz,ind]
          DivC[iP,jP,iz] = Div[iz,ind]
        end
      end
    end
    @views FRotCurl2Vec!(FCG[:,:,:,iF,uPos:vPos],Rot1C,Rot2C,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,iF,uPos:vPos],Grad1C,Grad2C,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,iF,ThPos],DivC,RhoCG,CG,Global,iF)
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])

    KE = TCacheCC1[Threads.threadid()]
    @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + wCCG*wCCG);

    @views FDiv3Vec!(FCG[:,:,:,iF,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

    if Global.Model.RefProfile
  #   Global.Pres=Pressure(ThCG,RhoCG,KE,Global)-Global.pBGrd;
  #   FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Global.Pres,CG,Global);
  #   FCG[:,:,:,:,uPos]=FCG[:,:,:,:,uPos]./RhoCG;
  #   FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
  #   if Global.Buoyancy
  #     RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
  #     FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
  #       Global.Grav*Global.JF[:,:,:,2:nz].*(RhoF-Global.RhoBGrdF)) ./ RhoF;
  #   end
    else
      Pres = TCacheCC2[Threads.threadid()]
      Pressure!(Pres,ThCG,RhoCG,KE,Global);
      @views FGrad3RhoVec!(FCG[:,:,:,iF,:],Pres,RhoCG,CG,Global,iF)
      if Global.Model.Buoyancy
        @inbounds for iz=1:nz-1  
          @inbounds for j=1:OP  
            @inbounds for i=1:OP  
              FCG[i,j,iz,iF,wPos] -= Grav*JF[i,j,iz+1,iF]
            end
          end  
        end
      end
    end
#   3-dim Curl and Grad of kinetic Energy
    @views FGrad3Vec!(FCG[:,:,:,iF,:],KE,CG,Global,iF)
    @views FCurlNon3Vec!(FCG[:,:,:,iF,:],v1CG,v2CG,wCG,wCCG,CG,Global,iF);

#   Divergence of Thermodynamic Variable
    if Global.Model.Thermo == "Energy"
      else
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,iF,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,iF,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    end  
  end  
  @inbounds for iF = 1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,iF,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,iF,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,iF,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,iF,ThPos] / CG.M[iz,ind]
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,iF,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  

  @inbounds for iG=1:CG.NumG
    if Global.Model.Damping
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
    if Global.Model.Source
      @views Source!(F[:,iG,:],U[:,iG,:],CG,Global);
    end
  end
end

function FcnNHCurlVec!(F,U,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

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
  Rot1 = Global.Cache.Rot1
  Rot2 = Global.Cache.Rot2
  Grad1 = Global.Cache.Grad1
  Grad2 = Global.Cache.Grad2
  Div = Global.Cache.Div
  DivTr = Global.Cache.DivTr
  FCG=Global.Cache.FCC
  Rot1CG=Global.Cache.Rot1C
  Rot2CG=Global.Cache.Rot2C
  Grad1CG=Global.Cache.Grad1C
  Grad2CG=Global.Cache.Grad2C
  DivCG=Global.Cache.DivC
  DivTrCG=Global.Cache.DivC
  RhoCG = Global.Cache.RhoCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  ThCG = Global.Cache.ThCG
  TrCG = Global.Cache.TrCG
  KE = Global.Cache.KE
  Pres = Global.Cache.Pres
  uStar = Global.Cache.uStar
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  DivTr .= 0.0
  F .= 0.0
  # Hyperdiffusion 
  @inbounds for iF = 1:NF
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

    #FTempW=FDivGrad2Vec(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,Global);
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
    for iT=1:NumTr
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

  @inbounds for iF = 1:NF
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
    @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + wCCG*wCCG);
#   Pressure
    @views Pres = Global.Cache.Pres[:,:,:,iF]  
    Pres1D = reshape(Pres,OP*OP*nz,1)
    Th1D = reshape(ThCG,OP*OP*nz,1)
    Rho1D = reshape(RhoCG,OP*OP*nz,1)
    Tr1D = reshape(TrCG,OP*OP*nz,NumTr)
#   @time Pressure!(reshape(v1CG,OP*OP*nz,1),reshape(ThCG,OP*OP*nz,1),reshape(RhoCG,OP*OP*nz,1),
#     reshape(TrCG,OP*OP*nz,NumTr),Global)
    Pressure!(Pres1D,Th1D,Rho1D,Tr1D,Global)
#   uStar
    @views uStarCoefficient!(uStar[:,:,iF],v1CG[:,:,1],v2CG[:,:,1],wCCG[:,:,1],CG,Global,iF)

#   Vertical Diffusion coefficient    
    KV = Global.Cache.DivC
    eddy_diffusivity_coefficient!(KV,v1CG,v2CG,wCCG,RhoCG,CG,Global,iF)

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
      @views FGrad3RhoVec!(FCG,Pres,RhoCG,CG,Global,iF)
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
        @views VerticalDiffusionSaclar!(FCG[:,:,:,iT+NumV],TrCG[:,:,:,iT],RhoCG,KV,CG,Global,iF)
        @views BoundaryFluxSaclar!(FCG[:,:,1,iT+NumV],TrCG[:,:,1,iT],Global.Cache.cTrS[:,:,iF,iT],CG,Global,iF)
      end  
    end

    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  

  if Global.Model.Damping
    @inbounds for iG=1:CG.NumG
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
  end
  if Global.Model.Source
    @inbounds for iG=1:CG.NumG
      @views Source!(F[:,iG,:],U[:,iG,:],CG,Global,iG)
    end
  end
  if Global.Model.Microphysics
    @inbounds for iG=1:CG.NumG
      @views SourceMicroPhysics(F[:,iG,:],U[:,iG,:],CG,Global,iG)
    end  
  end
end

function FcnNHCurlVecStatic!(F,U,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

(;  RhoPos,
    uPos,
    vPos,
    wPos,
    ThPos,
    NumV) = Global.Model

  Grav=Global.Phys.Grav    
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;
  JF = Global.Metric.JF
  Rot1 = Global.Cache.Rot1
  Rot2 = Global.Cache.Rot2
  Grad1 = Global.Cache.Grad1
  Grad2 = Global.Cache.Grad2
  Div = Global.Cache.Div
  FCG=Global.Cache.FCC
  Rot1CG=Global.Cache.Rot1C
  Rot2CG=Global.Cache.Rot2C
  Grad1CG=Global.Cache.Grad1C
  Grad2CG=Global.Cache.Grad2C
  DivCG=Global.Cache.DivC
  RhoCG = Global.Cache.RhoCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  ThCG = Global.Cache.ThCG
  KE = Global.Cache.KE
  Pres = Global.Cache.Pres
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  F .= 0.0
  # Hyperdiffusion 
  @inbounds for iF = 1:NF
    @views Ind = SVector{5*5}(reshape(CG.Glob[:,:,iF],5*5))
    @inbounds for iz=1:nz
      @views RhoCCG=SMatrix{5,5}(U[iz,Ind,RhoPos])
      @views v1CCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),uPos])
      @views v2CCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),vPos])
      @views ThCCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),ThPos])

      @views FRotCurl2VecDSS!(Rot1CG[:,:,iz],Rot2CG[:,:,iz],v1CCG,v2CCG,
        Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache)
      @views FGradDiv2VecDSS!(Grad1CG[:,:,iz],Grad2CG[:,:,iz],v1CCG,v2CCG,
        Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache,Val(5))

      #FTempW=FDivGrad2Vec(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,Global);
      @views FDivGrad2VecDSS!(DivCG[:,:,iz],ThCCG,RhoCCG,  
        Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache)

      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
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
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
        end
      end
    end
    @. FCG = 0.0
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF,Val(5))
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])

    @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + wCCG*wCCG);

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
      Pressure!(Pres,ThCG,RhoCG,KE,Global);
      @views FGrad3RhoVec!(FCG,Pres,RhoCG,CG,Global,iF)
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
    end  
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  

  @inbounds for iG=1:CG.NumG
    if Global.Model.Damping
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
    if Global.Model.Source
      @views Source!(F[:,iG,:],U[:,iG,:],CG,Global);
    end
  end
end

function FcnNHCurlVecTrStatic!(F,U,CG,Global)
  @unpack TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
  TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5 = Global.ThreadCache

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
  Rot1 = Global.Cache.Rot1
  Rot2 = Global.Cache.Rot2
  Grad1 = Global.Cache.Grad1
  Grad2 = Global.Cache.Grad2
  Div = Global.Cache.Div
  DivTr = Global.Cache.DivTr
  FCG=Global.Cache.FCC
  Rot1CG=Global.Cache.Rot1C
  Rot2CG=Global.Cache.Rot2C
  Grad1CG=Global.Cache.Grad1C
  Grad2CG=Global.Cache.Grad2C
  DivCG=Global.Cache.DivC
  RhoCG = Global.Cache.RhoCG
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  ThCG = Global.Cache.ThCG
  KE = Global.Cache.KE
  Pres = Global.Cache.Pres
  Rot1 .= 0.0
  Rot2 .= 0.0
  Grad1 .= 0.0
  Grad2 .= 0.0
  Div .= 0.0
  DivTr .= 0.0
  F .= 0.0
  # Hyperdiffusion 
  @time @inbounds for iF = 1:NF
    @views Ind = SVector{5*5}(reshape(CG.Glob[:,:,iF],5*5))
    @inbounds for iz=1:nz
      @views RhoCCG=SMatrix{5,5}(U[iz,Ind,RhoPos])
      @views v1CCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),uPos])
      @views v2CCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),vPos])
      @views ThCCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),ThPos])

      @views FRotCurl2VecDSS!(Rot1CG[:,:,iz],Rot2CG[:,:,iz],v1CCG,v2CCG,
        Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache)
      @time @views FGradDiv2VecDSS!(Grad1CG[:,:,iz],Grad2CG[:,:,iz],v1CCG,v2CCG,
        Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache,Val(5))

      @views FDivGrad2VecDSS!(DivCG[:,:,iz],ThCCG,RhoCCG,  
        Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache)

      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          Rot1[iz,ind] += Rot1CG[iP,jP,iz] / CG.M[iz,ind]
          Rot2[iz,ind] += Rot2CG[iP,jP,iz] / CG.M[iz,ind]
          Grad1[iz,ind] += Grad1CG[iP,jP,iz] / CG.M[iz,ind]
          Grad2[iz,ind] += Grad2CG[iP,jP,iz] / CG.M[iz,ind]
          Div[iz,ind] += DivCG[iP,jP,iz] / CG.M[iz,ind]
        end
      end
      for iT=1:NumTr
        @views ThCCG=SMatrix{5,5}(U[iz,reshape(CG.Glob[:,:,iF],5*5),iT+NumV])  
        @views FDivGrad2VecDSS!(DivCG[:,:,iz],ThCCG,RhoCCG,  
          Global.Metric.dXdxIC[:,:,iz,:,:,iF],Global.Metric.JC[:,:,iz,iF],CG,Global.ThreadCache)
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
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
          ThCG[iP,jP,iz] = U[iz,ind,ThPos]
          Rot1CG[iP,jP,iz] = Rot1[iz,ind]
          Rot2CG[iP,jP,iz] = Rot2[iz,ind]
          Grad1CG[iP,jP,iz] = Grad1[iz,ind]
          Grad2CG[iP,jP,iz] = Grad2[iz,ind]
          DivCG[iP,jP,iz] = Div[iz,ind]
        end
      end
    end
    @. FCG = 0.0
    @views FRotCurl2Vec!(FCG[:,:,:,uPos:vPos],Rot1CG,Rot2CG,CG,Global,iF)
    @views FGradDiv2Vec!(FCG[:,:,:,uPos:vPos],Grad1CG,Grad2CG,CG,Global,iF)
    @views FDivRhoGrad2Vec!(FCG[:,:,:,ThPos],DivCG,RhoCG,CG,Global,iF)
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
    @views @. wCCG = 0.5*(wCG[:,:,1:nz] + wCG[:,:,2:nz+1])

    @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + wCCG*wCCG);

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
      Pressure!(Pres,ThCG,RhoCG,KE,Global);
      @views FGrad3RhoVec!(FCG,Pres,RhoCG,CG,Global,iF)
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
    end  
    for iT = 1:NumTr
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            DivCG[iP,jP,iz] = DivTr[iz,ind,iT]
            ThCG[iP,jP,iz] = U[iz,ind,iT+NumV]
          end
        end
      end
      @views FDivRhoGrad2Vec!(FCG[:,:,:,iT+NumV],DivCG,RhoCG,CG,Global,iF)
      if Global.Model.Upwind
        @views FDiv3UpwindVec!(FCG[:,:,:,iT+NumV],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
      else
        @views FDiv3Vec!(FCG[:,:,:,iT+NumV],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
      end
    end  
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          F[iz,ind,RhoPos] += FCG[iP,jP,iz,RhoPos] / CG.M[iz,ind]
          F[iz,ind,uPos] += FCG[iP,jP,iz,uPos] / CG.M[iz,ind]
          F[iz,ind,vPos] += FCG[iP,jP,iz,vPos] / CG.M[iz,ind]
          F[iz,ind,ThPos] += FCG[iP,jP,iz,ThPos] / CG.M[iz,ind]
          for iT = 1:NumTr
            F[iz,ind,iT+NumV] += FCG[iP,jP,iz,iT+NumV] / CG.M[iz,ind]
          end  
        end
        @inbounds for iz=1:nz-1
          F[iz,ind,wPos] += FCG[iP,jP,iz,wPos] / CG.MW[iz,ind]
        end
      end  
    end
  end  

  @inbounds for iG=1:CG.NumG
    if Global.Model.Damping
      @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
    end
    if Global.Model.Source
      @views Source!(F[:,iG,:],U[:,iG,:],CG,Global);
    end
  end
end

