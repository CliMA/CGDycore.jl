function [F]=FcnNHCurlVec(U,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
RhoCG=reshape(U(reshape(CG.Glob,OP*OP*NF,1),:,RhoPos)...
  ,OP,OP,NF,nz);
v1CG=reshape(U(reshape(CG.Glob,OP*OP*NF,1),:,uPos)...
  ,OP,OP,NF,nz);
v2CG=reshape(U(reshape(CG.Glob,OP*OP*NF,1),:,vPos)...
  ,OP,OP,NF,nz);
wCG=zeros(OP,OP,NF,nz+1);
wCG(:,:,:,2:nz+1)=reshape(U(reshape(CG.Glob,OP*OP*NF,1),:,wPos)...
  ,OP,OP,NF,nz);
%wCG(:,:,:,1)=BoundaryW(v1CG,v2CG,CG,Param);
wCCG=0.5*(wCG(:,:,:,1:nz)+wCG(:,:,:,2:nz+1));
wCCG(:,:,:,1)=BoundaryW(v1CG,v2CG,CG,Param);
wCG(:,:,:,1)=2*wCCG(:,:,:,1)-wCG(:,:,:,2);
ThCG=reshape(U(reshape(CG.Glob,OP*OP*NF,1),:,ThPos)...
  ,OP,OP,NF,nz);

FCG=zeros(OP,OP,NF,nz,Param.NumV);
KE=0.5*(v1CG.*v1CG+v2CG.*v2CG+wCCG.*wCCG);
FCG(:,:,:,:,RhoPos)=-FDiv3Vec(RhoCG,v1CG,v2CG,wCG,CG,Param);

if Param.RefProfile
  Param.Pres=Pressure(ThCG,RhoCG,KE,Param)-Param.pBGrd;
  FCG(:,:,:,:,uPos:wPos)=-FGrad3Vec(Param.Pres,CG,Param);
  FCG(:,:,:,:,uPos)=FCG(:,:,:,:,uPos)./RhoCG;
  FCG(:,:,:,:,vPos)=FCG(:,:,:,:,vPos)./RhoCG;
  if Param.Buoyancy
    RhoF=0.5*(RhoCG(:,:,:,1:nz-1)+RhoCG(:,:,:,2:nz));
    FCG(:,:,:,1:nz-1,wPos)=(FCG(:,:,:,1:nz-1,wPos)...
      -Param.Grav*Param.JF(:,:,:,2:nz).*(RhoF-Param.RhoBGrdF))...
      ./RhoF;
  end
else
  Param.Pres=Pressure(ThCG,RhoCG,KE,Param);
  FCG(:,:,:,:,uPos:wPos)=-FGrad3Vec(Param.Pres,CG,Param);
  FCG(:,:,:,:,uPos)=FCG(:,:,:,:,uPos)./RhoCG;
  FCG(:,:,:,:,vPos)=FCG(:,:,:,:,vPos)./RhoCG;
  if Param.Buoyancy
    FCG(:,:,:,1:nz-1,wPos)=FCG(:,:,:,1:nz-1,wPos)...
      ./(0.5*(RhoCG(:,:,:,1:nz-1)+RhoCG(:,:,:,2:nz)))...
      -Param.Grav*Param.JF(:,:,:,2:nz);
  end
end
FCG(:,:,:,:,uPos:wPos)=FCG(:,:,:,:,uPos:wPos)-FGrad3Vec(KE,CG,Param);
FCG(:,:,:,:,uPos:wPos)=FCG(:,:,:,:,uPos:wPos)...
  +FCurlNon3Vec(v1CG,v2CG,wCG,wCCG,CG,Param);
if strcmp(Param.Thermo,'Energy')
else
  if Param.Upwind
    FCG(:,:,:,:,ThPos)=-FDiv3UpwindVec(ThCG,v1CG,v2CG,wCG,RhoCG,CG,Param);
  else
    FCG(:,:,:,:,ThPos)=-FDiv3Vec(ThCG,v1CG,v2CG,wCG,CG,Param);
  end
end
if Param.HyperVisc
  if strcmp(Param.Thermo,'Energy')
  else
    FCG(:,:,:,:,uPos:ThPos)=FCG(:,:,:,:,uPos:ThPos)...
      +HyperDiffusionVec(v1CG,v2CG,wCG,ThCG,RhoCG,CG,Param);
  end
end



F=zeros(size(U));
for iM=1:size(CG.FaceGlob,2)
  F(reshape(CG.Glob(:,CG.FaceGlob(iM).Ind,:)...
    ,OP*OP*size(CG.FaceGlob(iM).Ind,1),1),:,:)...
    =F(reshape(CG.Glob(:,CG.FaceGlob(iM).Ind,:)...
    ,OP*OP*size(CG.FaceGlob(iM).Ind,1),1),:,:)...
    +reshape(FCG(:,:,CG.FaceGlob(iM).Ind,:)...
    ,OP*OP*size(CG.FaceGlob(iM).Ind,1),nz,Param.NumV);
end
F(:,:,RhoPos)=F(:,:,RhoPos)./CG.M;
F(:,:,uPos)=F(:,:,uPos)./CG.M;
F(:,:,vPos)=F(:,:,vPos)./CG.M;
F(:,1:nz-1,wPos)=F(:,1:nz-1,wPos)./CG.MW;
F(:,:,ThPos)=F(:,:,ThPos)./CG.M;

if Param.Damping
  F(:,1:nz-1,wPos)=F(:,1:nz-1,wPos)+Damping(U(:,:,wPos),Param);
end
if Param.Source
  F=F+Source(U,CG,Param);
end
end

function F=Damping(W,Param)
F=zeros(size(W,1),size(W,2)-2);
for iz=Param.Grid.nz-1:-1:1
  zLoc=Param.Grid.z(iz+1);
  if zLoc>=Param.H-Param.StrideDamp
    Damp = Param.Relax*...
      sin(0.5*pi*(1.0 - (Param.H - zLoc)/Param.StrideDamp))^2;
    F(:,iz)=-Damp*W(:,iz);
  else
    break
  end
end
end

