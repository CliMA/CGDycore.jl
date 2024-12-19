function [fig] = PlotDG(y,Grid,Trans,fig,Param)
OrdPoly=Param.OrdPolyX;
figure(fig)
clf(fig)
hold on;
for iF=1:Grid.NumFaces
  X=zeros(size(Grid.Faces(iF).N,2),3);
  cLoc=reshape(y(1,1:OrdPoly+1,1:OrdPoly+1,iF),OrdPoly+1,OrdPoly+1);
  if strcmp(Grid.Type,'Quad')
    dd=2/OrdPoly;
    ksi0=-1;
    for iRef=1:OrdPoly
      ksi1=ksi0+dd;
      eta0=-1;
      for jRef=1:OrdPoly
        eta1=eta0+dd;
        c=0;
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            c=c+Lagrange(0.5*(ksi0+ksi1),Param.xw,i)*Lagrange(0.5*(eta0+eta1),Param.xw,j)...
            *cLoc(i,j);
          end
        end
        [X(1,:),det,DF]=Trans(ksi0,eta0,Grid.Faces(iF),Grid);
        [X(2,:),det,DF]=Trans(ksi1,eta0,Grid.Faces(iF),Grid);
        [X(3,:),det,DF]=Trans(ksi1,eta1,Grid.Faces(iF),Grid);
        [X(4,:),det,DF]=Trans(ksi0,eta1,Grid.Faces(iF),Grid);
        fill3(X(:,1),X(:,2),X(:,3),c);
        eta0=eta1;
      end
      ksi0=ksi1;
    end
  else
  end
end
fig=fig+1;
hold off
end

