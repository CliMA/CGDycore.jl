function [fig] = Plot2(c,Trans,Param)
OrdPoly=Param.OrdPoly;
Subs=Param.OrdPoly;
fig=Param.fig;
wQX=Param.w;
xQX=Param.xw;
wQZ=Param.w;
xQZ=Param.xw;
figure(fig);
hold on
dS=2/Subs;
for iF=1:size(Param.Grid.Faces,2)
  ksi=-1;
  for is=1:Subs
    eta=-1;
    for js=1:Subs
      cP=0;
      for jx=1:OrdPoly+1
        for jz=1:OrdPoly+1
          cP=cP+Lagrange(ksi+0.5*dS,xQX,jx)*Lagrange(eta+0.5*dS,xQZ,jz)...
            *c(1,jx,jz,iF);
        end
      end
      [X1]=Trans(ksi,eta,Param.Grid.Faces(iF),Param.Grid);
      [X2]=Trans(ksi+dS,eta,Param.Grid.Faces(iF),Param.Grid);
      [X3]=Trans(ksi+dS,eta+dS,Param.Grid.Faces(iF),Param.Grid);
      [X4]=Trans(ksi,eta+dS,Param.Grid.Faces(iF),Param.Grid);
      
      if strcmp(Param.GridType,'Radial')
        r1=sqrt(X1(1,1,1)^2+X1(1,1,2)^2);
        X1=(Param.RadIP+(r1-Param.RadI)/(Param.RadO-Param.RadI))/r1*X1;
        r2=sqrt(X2(1,1,1)^2+X2(1,1,2)^2);
        X2=(Param.RadIP+(r2-Param.RadI)/(Param.RadO-Param.RadI))/r2*X2;
        r3=sqrt(X3(1,1,1)^2+X3(1,1,2)^2);
        X3=(Param.RadIP+(r3-Param.RadI)/(Param.RadO-Param.RadI))/r3*X3;
        r4=sqrt(X4(1,1,1)^2+X4(1,1,2)^2);
        X4=(Param.RadIP+(r4-Param.RadI)/(Param.RadO-Param.RadI))/r4*X4;
      end
      x=[X1(1) X2(1) X3(1) X4(1) X1(1)];
      y=[X1(2) X2(2) X3(2) X4(2) X1(2)];
      
      fill(x,y,cP);
      eta=eta+dS;
    end
    ksi=ksi+dS;
  end
end
hold off
fig=fig+1;
end

