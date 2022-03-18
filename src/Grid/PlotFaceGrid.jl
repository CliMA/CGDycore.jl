function []=PlotFaceGrid(Grid,Subs,Trans,fig)
figure(fig)
clf(fig)
hold on;
dS=2/Subs;
ksi=zeros(2,1);
for iF=1:size(Grid.Faces,2)
  ksi(1)=-1;
  for is=1:Subs
    ksi(2)=-1;
    for js=1:Subs
      [J,D,X1]=Trans(ksi(1),ksi(2),Grid.Faces(iF),Grid);
      [J,D,X2]=Trans(ksi(1)+dS,ksi(2),Grid.Faces(iF),Grid);
      [J,D,X3]=Trans(ksi(1)+dS,ksi(2)+dS,Grid.Faces(iF),Grid);
      [J,D,X4]=Trans(ksi(1),ksi(2)+dS,Grid.Faces(iF),Grid);
      
      x=[X1(1) X2(1) X3(1) X4(1) X1(1)];
      y=[X1(2) X2(2) X3(2) X4(2) X1(2)];
      if size(X1,3)==2
        fill(x,y,1);
      else
        z=[X1(3) X2(3) X3(3) X4(3) X1(3)];
        fill3(x,y,z,1);
      end
      ksi(2)=ksi(2)+dS;
    end
    ksi(1)=ksi(1)+dS;
  end
end
hold off
end
