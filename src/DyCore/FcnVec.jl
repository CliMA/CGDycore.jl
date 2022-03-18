function [F]=FcnVec(U,CG,Param)
nz=Param.Grid.nz;
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
W=[zeros(size(U,1),1)...
   ,U(:,:,wPos)];
W(:,1)=BoundaryW(U(:,:,uPos),U(:,:,vPos),CG,Param);


WC=0.5*(W(:,1:nz)+W(:,2:nz+1));
KE=0.5*(U(:,:,Param.uPos).*U(:,:,Param.uPos)+U(:,:,Param.vPos).*U(:,:,Param.vPos)+WC.*WC);
FRho=-FDiv3Vec(U(:,:,RhoPos),U(:,:,uPos:vPos),W,CG,Param);
Param.Pres=Pressure(U(:,:,ThPos),U(:,:,RhoPos),KE,Param);
[FUG,FWG]=FGrad3Vec(Param.Pres,CG,Param);
FUG(:,:,1)=FUG(:,:,1)./U(:,:,RhoPos);
FUG(:,:,2)=FUG(:,:,2)./U(:,:,RhoPos);
if Param.Buoyancy
  FWG(:,:)=FWG(:,:)./(0.5*(U(:,1:nz-1,RhoPos)...
    +U(:,2:nz,RhoPos)))+Param.Grav;
else
  FWG(:,:)=FWG(:,:)./(0.5*(U(:,1:nz-1,RhoPos)+U(:,2:nz,RhoPos)));
end
[FUKE,FWKE]=FGrad3Vec(KE,CG,Param);

[FU2,FW2]=FCurlNon3Vec(U(:,:,uPos:vPos),W,CG,Param);
FU=-FUG-FUKE+FU2;
FW=-FWG-FWKE+FW2;
if strcmp(Param.Thermo,'Energy')
  FTh=-FDiv3Vec(U(:,:,ThPos)+Param.Pres,U(:,:,uPos:vPos),W,CG,Param);
else
  FTh=-FDiv3Vec(U(:,:,ThPos),U(:,:,uPos:vPos),W,CG,Param);
end
if Param.HyperVisc
  if strcmp(Param.Thermo,'Energy')
    [FUD,FWD,FTD]=HyperDiffusionVec(U(:,:,uPos:vPos),W...
      ,U(:,:,ThPos)+Param.Pres...
      ,U(:,:,RhoPos),CG,Param);
  else
    [FUD,FWD,FTD]=HyperDiffusionVec(U(:,:,uPos:vPos),W,U(:,:,ThPos)...
      ,U(:,:,RhoPos),CG,Param);
  end
  FU=FU+FUD;
  %FW=FW+FWD;
  FTh=FTh+FTD;
end
if Param.Damping
  FW=FW+Damping(W,Param);
end
if Param.Source
  [FUD,FTD]=Source(U,CG,Param);
  FU=FU+FUD;
  FTh=FTh+FTD;
end
  
F=zeros(size(U));

F(:,:,RhoPos)=FRho;
F(:,:,uPos:vPos)=FU;
F(:,1:nz-1,wPos)=FW;
F(:,:,ThPos)=FTh;
end

function F=Damping(W,Param)
F=zeros(size(W,1),size(W,2)-2);
for iz=Param.Grid.nz-1:-1:1
  zLoc=Param.Grid.z(iz+1);
  if zLoc>=Param.H-Param.StrideDamp
    Damp = Param.Relax*...
      sin(0.5*pi*(1.0 - (Param.H - zLoc)/Param.StrideDamp))^2;
    F(:,iz)=-Damp*W(:,iz+1);
  else
    break
  end
end
end

