function [cP]=DG2Project(c,Param)
OrdPolyX=Param.OrdPolyX;
NumFaces = size(Param.X,4);
cP=zeros(1,OrdPolyX+1,OrdPolyX+1,NumFaces);
for iF=1:NumFaces
  for i=1:OrdPolyX+1
    for j=1:OrdPolyX+1
      cP(:,i,j,iF)=c(Param.X(i,j,:,iF),Param);
    end
  end
end
end

