function FcnNHCurlVec!(F,U,CG,Global)
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
FCG = Global.Cache.FCG
FCG .= 0
KE = Global.Cache.KE
Pres = Global.Cache.Pres
RhoCG = Global.Cache.RhoCG
v1CG = Global.Cache.v1CG
v2CG = Global.Cache.v2CG
wCG = Global.Cache.wCG
wCCG = Global.Cache.wCCG
ThCG = Global.Cache.ThCG
JF = Global.Metric.JF
Rot1CG = Global.Cache.Rot1CG
Rot2CG = Global.Cache.Rot2CG
Grad1CG = Global.Cache.Grad1CG
Grad2CG = Global.Cache.Grad2CG
DivCG = Global.Cache.DivCG
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
@inbounds for iF=1:NF
  iG=0
  @inbounds for jP=1:OP
    @inbounds for iP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        RhoCG[iP,jP,iz] = U[iz,ind,RhoPos]
        v1CG[iP,jP,iz] = U[iz,ind,uPos]
        v2CG[iP,jP,iz] = U[iz,ind,vPos]
        ThCG[iP,jP,iz] = U[iz,ind,ThPos]
      end
    end
  end

  FRotCurl2VecDSS!(Rot1CG,Rot2CG,v1CG,v2CG,CG,Global,iF)
  FGradDiv2VecDSS!(Grad1CG,Grad2CG,v1CG,v2CG,CG,Global,iF)

# FTempW=FDivGrad2Vec(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,Global);
  FDivGrad2VecDSS!(DivCG,ThCG,RhoCG,CG,Global,iF)  
  iG=0
  @inbounds for jP=1:OP
    @inbounds for iP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        Rot1[iz,ind] += Rot1CG[iP,jP,iz]
        Rot2[iz,ind] += Rot2CG[iP,jP,iz]
        Grad1[iz,ind] += Grad1CG[iP,jP,iz]
        Grad2[iz,ind] += Grad2CG[iP,jP,iz]
        Div[iz,ind] += DivCG[iP,jP,iz]
      end
    end
  end
end

@. Rot1 = Rot1 / CG.M
@. Rot2 = Rot2 / CG.M
@. Grad1 = Grad1 / CG.M
@. Grad2 = Grad2 / CG.M
@. Div = Div / CG.M

@inbounds for iF=1:NF
  FCG .= 0.0  
  iG=0
  @inbounds for jP=1:OP
    @inbounds for iP=1:OP
      iG=iG+1
      ind=CG.Glob[iG,iF]
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
  @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + wCCG*wCCG);

  @views FDiv3Vec!(FCG[:,:,:,RhoPos],RhoCG,v1CG,v2CG,wCG,CG,Global,iF);

  if Global.Model.RefProfile
  # Global.Pres=Pressure(ThCG,RhoCG,KE,Global)-Global.pBGrd;
  # FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Global.Pres,CG,Global);
  # FCG[:,:,:,:,uPos]=FCG[:,:,:,:,uPos]./RhoCG;
  # FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
  # if Global.Buoyancy
  #   RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
  #   FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
  #     Global.Grav*Global.JF[:,:,:,2:nz].*(RhoF-Global.RhoBGrdF)) ./ RhoF;
  # end
  else
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
# 3-dim Curl and Grad of kinetic Energy
  FGrad3Vec!(FCG,KE,CG,Global,iF)
  FCurlNon3Vec!(FCG,v1CG,v2CG,wCG,wCCG,CG,Global,iF);

# Divergence of Thermodynamic Variable
  if Global.Model.Thermo == "Energy"
  else
    if Global.Model.Upwind
      @views FDiv3UpwindVec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global,iF);
    else
      @views FDiv3Vec!(FCG[:,:,:,ThPos],ThCG,v1CG,v2CG,wCG,CG,Global,iF);
    end
  end  
  iG=0
  @inbounds for jP=1:OP
    @inbounds for iP=1:OP
      iG=iG+1  
      ind=CG.Glob[iG,iF]
      @inbounds for iz=1:nz
        @inbounds for iV=1:NumV  
          F[iz,ind,iV] += FCG[iP,jP,iz,iV]  
        end
      end
    end  
  end
end  
@views @. F[:,:,RhoPos] = F[:,:,RhoPos] / CG.M;
@views @. F[:,:,uPos] = F[:,:,uPos] / CG.M;
@views @. F[:,:,vPos] = F[:,:,vPos] / CG.M;
@views @. F[1:nz-1,:,wPos] = F[1:nz-1,:,wPos] / CG.MW;
@views @. F[:,:,ThPos] = F[:,:,ThPos] / CG.M;

for iG=1:CG.NumG
  if Global.Model.Damping
    @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
  end
end
if Global.Model.Source
  Source!(F,U,CG,Global);
end
end

function FcnNHSourceVec!(F,U,CG,Global)
ThPos=Global.ThPos;
F .= 0
if Global.Source
  Source!(F,U,CG,Global);
end
end

function FcnNHCurlVec(U,CG,Global)
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;
RhoPos=Global.RhoPos;
uPos=Global.uPos;
vPos=Global.vPos;
wPos=Global.wPos;
ThPos=Global.ThPos;
RhoCG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,RhoPos] ,OP,OP,NF,nz);
v1CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,uPos] ,OP,OP,NF,nz);
v2CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,vPos] ,OP,OP,NF,nz);
wCG=zeros(OP,OP,NF,nz+1);
wCG[:,:,:,2:nz+1]=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,wPos] ,OP,OP,NF,nz);
wCCG=0.5*(wCG[:,:,:,1:nz]+wCG[:,:,:,2:nz+1]);
wCCG[:,:,:,1]=BoundaryW(v1CG,v2CG,CG,Global);
wCG[:,:,:,1]=2*wCCG[:,:,:,1]-wCG[:,:,:,2];
ThCG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,ThPos]
  ,OP,OP,NF,nz);

FCG=zeros(OP,OP,NF,nz,Global.NumV);
KE=0.5*(v1CG.*v1CG+v2CG.*v2CG+wCCG.*wCCG);

FDiv3Vec!(view(FCG,:,:,:,:,RhoPos),RhoCG,v1CG,v2CG,wCG,CG,Global);
if Global.RefProfile
# Global.Pres=Pressure(ThCG,RhoCG,KE,Global)-Global.pBGrd;
# FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Global.Pres,CG,Global);
# FCG[:,:,:,:,uPos]=FCG[:,:,:,:,uPos]./RhoCG;
# FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
# if Global.Buoyancy
#   RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
#   FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
#     Global.Grav*Global.JF[:,:,:,2:nz].*(RhoF-Global.RhoBGrdF)) ./ RhoF;
# end
else
  Global.Pres=Pressure(ThCG,RhoCG,KE,Global);
  FGrad3RhoVec!(FCG,Global.Pres,RhoCG,CG,Global)
  if Global.Buoyancy
    @views FCG[:,:,:,1:nz-1,wPos] .-= Global.Grav*Global.JF[:,:,:,2:nz];
  end
end
FGrad3Vec!(FCG,KE,CG,Global)
FCG[:,:,:,:,uPos:wPos]=FCG[:,:,:,:,uPos:wPos]+
  FCurlNon3Vec(v1CG,v2CG,wCG,wCCG,CG,Global);
if strcmp(Global.Thermo,"Energy")
else
  if Global.Upwind
    FCG[:,:,:,:,ThPos]=-FDiv3UpwindVec(ThCG,v1CG,v2CG,wCG,RhoCG,CG,Global);
  else
    FDiv3Vec!(view(FCG,:,:,:,:,ThPos),ThCG,v1CG,v2CG,wCG,CG,Global);
  end
end
if Global.HyperVisc
  if strcmp(Global.Thermo,"Energy")
  else
#   FCG[:,:,:,:,uPos:ThPos]=FCG[:,:,:,:,uPos:ThPos]+
#   HyperDiffusionVec(v1CG,v2CG,wCG,ThCG,RhoCG,CG,Global);
    HyperDiffusionVec!(FCG,v1CG,v2CG,wCG,ThCG,RhoCG,CG,Global)  
  end
end

F=zeros(size(U));
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  mat = reshape(FCG[:,:,i,:,:] ,OP*OP*size(i,1),nz,Global.NumV)
  F[arr,:,:] = F[arr,:,:] .+ mat;
end
F[:,:,RhoPos]=F[:,:,RhoPos]./CG.M;
F[:,:,uPos]=F[:,:,uPos]./CG.M;
F[:,:,vPos]=F[:,:,vPos]./CG.M;
F[:,1:nz-1,wPos]=F[:,1:nz-1,wPos]./CG.MW;
F[:,:,ThPos]=F[:,:,ThPos]./CG.M;


for iG=1:CG.NumG
  if Global.Model.Damping
    @views Damping!(F[:,iG,wPos],U[:,iG,wPos],Global)  
  end
end
if Global.Source
  Source!(F,U,CG,Global);
end
return F
end

function Damping!(F,W,Global)
  @inbounds for iz=Global.Grid.nz-1:-1:1
    zLoc=Global.Grid.z[iz];
    if zLoc>=Global.Grid.H-Global.Model.StrideDamp
      Damp = Global.Model.Relax*
        sin(0.5*pi*(1.0 - (Global.Grid.H - zLoc)/Global.Model.StrideDamp))^2;
      F[iz] -=Damp*W[iz];
    else
      break
    end
  end
end

