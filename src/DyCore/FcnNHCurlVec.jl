function FcnNHCurlVec!(F,U,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
FCG = Param.FCG
FCG .= 0
KE = Param.KE
Pres = Param.Pres
RhoCG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,RhoPos]
  ,OP,OP,NF,nz);
v1CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,uPos]
  ,OP,OP,NF,nz);
v2CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,vPos]
  ,OP,OP,NF,nz);
wCG=zeros(OP,OP,NF,nz+1);
wCG[:,:,:,2:nz+1]=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,wPos]
  ,OP,OP,NF,nz);
wCCG=0.5*(wCG[:,:,:,1:nz]+wCG[:,:,:,2:nz+1]);
wCCG[:,:,:,1]=BoundaryW(v1CG,v2CG,CG,Param);
wCG[:,:,:,1]=2*wCCG[:,:,:,1]-wCG[:,:,:,2];
ThCG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,ThPos]
  ,OP,OP,NF,nz);

KE .= 0.5 .* (v1CG.*v1CG+v2CG.*v2CG+wCCG.*wCCG);

FDiv3Vec!(view(FCG,:,:,:,:,RhoPos),RhoCG,v1CG,v2CG,wCG,CG,Param);
if Param.RefProfile
# Param.Pres=Pressure(ThCG,RhoCG,KE,Param)-Param.pBGrd;
# FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Param.Pres,CG,Param);
# FCG[:,:,:,:,uPos]=FCG[:,:,:,:,uPos]./RhoCG;
# FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
# if Param.Buoyancy
#   RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
#   FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
#     Param.Grav*Param.JF[:,:,:,2:nz].*(RhoF-Param.RhoBGrdF)) ./ RhoF;
# end
else
  Pres .= Pressure(ThCG,RhoCG,KE,Param);
  FGrad3RhoVec!(FCG,Pres,RhoCG,CG,Param)
  if Param.Buoyancy
    @views FCG[:,:,:,1:nz-1,wPos] .-= Param.Grav*Param.JF[:,:,:,2:nz];
  end
end
FGrad3Vec!(FCG,KE,CG,Param)
FCG[:,:,:,:,uPos:wPos]=FCG[:,:,:,:,uPos:wPos]+
  FCurlNon3Vec(v1CG,v2CG,wCG,wCCG,CG,Param);
if strcmp(Param.Thermo,"Energy")
else
  if Param.Upwind
    FCG[:,:,:,:,ThPos]=-FDiv3UpwindVec(ThCG,v1CG,v2CG,wCG,RhoCG,CG,Param);
  else
    FDiv3Vec!(view(FCG,:,:,:,:,ThPos),ThCG,v1CG,v2CG,wCG,CG,Param);
  end
end
if Param.HyperVisc
  if strcmp(Param.Thermo,"Energy")
  else
#   FCG[:,:,:,:,uPos:ThPos]=FCG[:,:,:,:,uPos:ThPos]+
#   HyperDiffusionVec(v1CG,v2CG,wCG,ThCG,RhoCG,CG,Param);
    HyperDiffusionVec!(FCG,v1CG,v2CG,wCG,ThCG,RhoCG,CG,Param)  
  end
end

F .= 0
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  mat = reshape(FCG[:,:,i,:,:] ,OP*OP*size(i,1),nz,Param.NumV)
  @views F[arr,:,:] = F[arr,:,:] .+ mat;
end
@views F[:,:,RhoPos]=F[:,:,RhoPos]./CG.M;
@views F[:,:,uPos]=F[:,:,uPos]./CG.M;
@views F[:,:,vPos]=F[:,:,vPos]./CG.M;
@views F[:,1:nz-1,wPos]=F[:,1:nz-1,wPos]./CG.MW;
@views F[:,:,ThPos]=F[:,:,ThPos]./CG.M;

if Param.Damping
  @views F[:,1:nz-1,wPos]=F[:,1:nz-1,wPos]+Damping(U[:,:,wPos],Param);
end
if Param.Source
  Source!(F,U,CG,Param);
end
end

function FcnNHCurlVec(U,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
RhoCG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,RhoPos]
  ,OP,OP,NF,nz);
v1CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,uPos]
  ,OP,OP,NF,nz);
v2CG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,vPos]
  ,OP,OP,NF,nz);
wCG=zeros(OP,OP,NF,nz+1);
wCG[:,:,:,2:nz+1]=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,wPos]
  ,OP,OP,NF,nz);
wCCG=0.5*(wCG[:,:,:,1:nz]+wCG[:,:,:,2:nz+1]);
wCCG[:,:,:,1]=BoundaryW(v1CG,v2CG,CG,Param);
wCG[:,:,:,1]=2*wCCG[:,:,:,1]-wCG[:,:,:,2];
ThCG=reshape(U[reshape(CG.Glob,OP*OP*NF,1),:,ThPos]
  ,OP,OP,NF,nz);

FCG=zeros(OP,OP,NF,nz,Param.NumV);
KE=0.5*(v1CG.*v1CG+v2CG.*v2CG+wCCG.*wCCG);

FDiv3Vec!(view(FCG,:,:,:,:,RhoPos),RhoCG,v1CG,v2CG,wCG,CG,Param);
if Param.RefProfile
# Param.Pres=Pressure(ThCG,RhoCG,KE,Param)-Param.pBGrd;
# FCG[:,:,:,:,uPos:wPos]=-FGrad3Vec(Param.Pres,CG,Param);
# FCG[:,:,:,:,uPos]=FCG[:,:,:,:,uPos]./RhoCG;
# FCG[:,:,:,:,vPos]=FCG[:,:,:,:,vPos]./RhoCG;
# if Param.Buoyancy
#   RhoF=0.5*(RhoCG[:,:,:,1:nz-1]+RhoCG[:,:,:,2:nz]);
#   FCG[:,:,:,1:nz-1,wPos]=(FCG[:,:,:,1:nz-1,wPos]-
#     Param.Grav*Param.JF[:,:,:,2:nz].*(RhoF-Param.RhoBGrdF)) ./ RhoF;
# end
else
  Param.Pres=Pressure(ThCG,RhoCG,KE,Param);
  FGrad3RhoVec!(FCG,Param.Pres,RhoCG,CG,Param)
  if Param.Buoyancy
    @views FCG[:,:,:,1:nz-1,wPos] .-= Param.Grav*Param.JF[:,:,:,2:nz];
  end
end
FGrad3Vec!(FCG,KE,CG,Param)
FCG[:,:,:,:,uPos:wPos]=FCG[:,:,:,:,uPos:wPos]+
  FCurlNon3Vec(v1CG,v2CG,wCG,wCCG,CG,Param);
if strcmp(Param.Thermo,"Energy")
else
  if Param.Upwind
    FCG[:,:,:,:,ThPos]=-FDiv3UpwindVec(ThCG,v1CG,v2CG,wCG,RhoCG,CG,Param);
  else
    FDiv3Vec!(view(FCG,:,:,:,:,ThPos),ThCG,v1CG,v2CG,wCG,CG,Param);
  end
end
if Param.HyperVisc
  if strcmp(Param.Thermo,"Energy")
  else
#   FCG[:,:,:,:,uPos:ThPos]=FCG[:,:,:,:,uPos:ThPos]+
#   HyperDiffusionVec(v1CG,v2CG,wCG,ThCG,RhoCG,CG,Param);
    HyperDiffusionVec!(FCG,v1CG,v2CG,wCG,ThCG,RhoCG,CG,Param)  
  end
end

F=zeros(size(U));
for iM=1:size(CG.FaceGlob,1)
  i = CG.FaceGlob[iM].Ind
  arr = reshape(CG.Glob[:,i,:],OP*OP*size(i,1))
  mat = reshape(FCG[:,:,i,:,:] ,OP*OP*size(i,1),nz,Param.NumV)
  F[arr,:,:] = F[arr,:,:] .+ mat;
end
F[:,:,RhoPos]=F[:,:,RhoPos]./CG.M;
F[:,:,uPos]=F[:,:,uPos]./CG.M;
F[:,:,vPos]=F[:,:,vPos]./CG.M;
F[:,1:nz-1,wPos]=F[:,1:nz-1,wPos]./CG.MW;
F[:,:,ThPos]=F[:,:,ThPos]./CG.M;

if Param.Damping
  F[:,1:nz-1,wPos]=F[:,1:nz-1,wPos]+Damping(U[:,:,wPos],Param);
end
if Param.Source
  Source!(F,U,CG,Param);
end
return F
end

function Damping(W,Param)
F=zeros(size(W,1),size(W,2)-1);
for iz=Param.Grid.nz-1:-1:1
  zLoc=Param.Grid.z[iz];
  if zLoc>=Param.H-Param.StrideDamp
    Damp = Param.Relax*
      sin(0.5*pi*(1.0 - (Param.H - zLoc)/Param.StrideDamp))^2;
    F[:,iz]=-Damp*W[:,iz];
  else
    break
  end
end
return F
end

