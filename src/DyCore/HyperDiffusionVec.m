function [FUCG]=HyperDiffusionVec(v1CG,v2CG,wCG,ThCG,RhoCG,CG,Param)
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
FUCG=zeros(OP,OP,NF,nz,4);
[FRTemp1,FRTemp2]=FRotCurl2VecDSS(v1CG,v2CG,CG,Param);
[FGTemp1,FGTemp2]=FGradDiv2VecDSS(v1CG,v2CG,CG,Param);
FRTemp1=FRTemp1-FGTemp1;
FRTemp2=FRTemp2-FGTemp2;
FUCG(:,:,:,:,1:2)=-Param.HyperDCurl*FRotCurl2Vec(FRTemp1,FRTemp2,CG,Param);
FUCG(:,:,:,:,1:2)=FUCG(:,:,:,:,1:2)...
  +Param.HyperDGrad*FGradDiv2Vec(FRTemp1,FRTemp2,CG,Param);

% FTempW=FDivGrad2Vec(0.5*(W(:,1:nz)+W(:,2:nz+1)),CG,Param);
% FTempW=-Param.HyperDDiv*FDivGrad2Vec(FTempW,CG,Param);
%FW(:,2:nz)=0.5*(FTempW(:,1:nz-1)+FTempW(:,2:nz));
% 
% FTempT=FDivGrad2VecDSS(ThCG./RhoCG,CG,Param);
% FUCG(:,:,:,:,4)=-Param.HyperDDiv*FDivRhoGrad2Vec(FTempT,RhoCG,CG,Param);
if Param.HyperDDiv>0
  FTempT=FDivGrad2VecDSS(ThCG./RhoCG,CG,Param);
  FUCG(:,:,:,:,4)=-Param.HyperDDiv*FDivRhoGrad2Vec(FTempT,RhoCG,CG,Param);
end
end
