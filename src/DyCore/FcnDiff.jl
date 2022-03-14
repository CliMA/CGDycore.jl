function [F]=FcnDiff(U,CG,Param)
nz=Param.Grid.nz;
RhoPos=Param.RhoPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
ThPos=Param.ThPos;
W=[zeros(size(U,1),1)...
   ,U(:,:,wPos)];

if Param.HyperD>0
  if strcmp(Param.Thermo,'Energy')
    [FU,FW,FTh]=DiffusionVec(U(:,:,uPos:vPos),W...
      ,U(:,:,ThPos)+Param.Pres...
      ,U(:,:,RhoPos),CG,Param);
  else
    [FU,FW,FTh]=DiffusionVec(U(:,:,uPos:vPos),W,U(:,:,ThPos)...
      ,U(:,:,RhoPos),CG,Param);
  end
end
 
F=zeros(size(U));
F(:,:,uPos:vPos)=FU;
F(:,:,wPos)=FW(:,2:nz+1);
F(:,:,ThPos)=FTh;
end

