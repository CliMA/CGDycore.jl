function [M,MW]=MassCG(CG,Param)
OrdPoly=CG.OrdPoly;
nz=Param.Grid.nz;
M=zeros(CG.NumG,nz);
for iF=1:Param.Grid.NumFaces
  M(CG.Faces(iF).Glob,:)=M(CG.Faces(iF).Glob,:)...
    +reshape(0.5*(Param.J(:,:,1,iF,:)+Param.J(:,:,2,iF,:))...
    ,(OrdPoly+1)*(OrdPoly+1),nz);
end
MW=0.5*(M(:,1:end-1)+M(:,2:end));
end
