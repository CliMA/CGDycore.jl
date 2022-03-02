function [fig] = PlotDG(c,CG,Trans,Topo,Param,fig,Slice)
OrdPoly=CG.OrdPoly;
figure(fig)
clf(fig)
hold on;
if strcmp(Param.Grid.Form,'Sphere') && Param.Flat
  dTol=2*pi/max(Param.nPanel-1,1);
end
switch Slice.Type
  case 'XY'
    level=Slice.iz;
    for iF=1:Param.Grid.NumFaces
      X=zeros(size(Param.Grid.Faces(iF).N,2),3);
      cLoc=reshape(c(:,:,iF),OrdPoly+1,OrdPoly+1);
      if strcmp(Param.Grid.Type,'Quad')
        dd=2/OrdPoly;
        ksi0=-1;
        for iRef=1:OrdPoly
          ksi1=ksi0+dd;
          eta0=-1;
          for jRef=1:OrdPoly
            eta1=eta0+dd;
            cc=0;
            for j=1:OrdPoly+1
              for i=1:OrdPoly+1
                cc=cc+Lagrange(0.5*(ksi0+ksi1),CG.xw,i)*Lagrange(0.5*(eta0+eta1),CG.xw,j)...
                  *cLoc(i,j);
              end
            end
            
            [X(1,:)]=Trans(ksi0,eta0,Param.Grid.z(level)...
              ,Param.Grid.Faces(iF),Topo,Param);
            [X(2,:)]=Trans(ksi1,eta0,Param.Grid.z(level)...
              ,Param.Grid.Faces(iF),Topo,Param);
            [X(3,:)]=Trans(ksi1,eta1,Param.Grid.z(level)...
              ,Param.Grid.Faces(iF),Topo,Param);
            [X(4,:)]=Trans(ksi0,eta1,Param.Grid.z(level)...
              ,Param.Grid.Faces(iF),Topo,Param);
            
            if strcmp(Param.Grid.Form,'Sphere') && Param.Flat
              [lam(1),theta(1)]=cart2sphere(X(1,1),X(1,2),X(1,3));
              [lam(2),theta(2)]=cart2sphere(X(2,1),X(2,2),X(2,3));
              [lam(3),theta(3)]=cart2sphere(X(3,1),X(3,2),X(3,3));
              [lam(4),theta(4)]=cart2sphere(X(4,1),X(4,2),X(4,3));
              lammin=min(lam);
              lammax=max(lam);
              if abs(lammin-lammax)>2*pi-dTol
                for i=1:4
                  if lam(i)<pi
                    lam(i)=lam(i)+2*pi;
                    if lam(i)>3*pi
                      lam(i)=lam(i)-2*pi;
                    end
                  end
                end
              end
              fill(lam,theta,real(cc));
            else
              fill3(X(:,1),X(:,2),X(:,3),real(cc));
            end
            
            eta0=eta1;
          end
          ksi0=ksi1;
        end
      else
      end
    end
  case 'XZ'
    y=Slice.y;
    for iF=1:Param.Grid.NumFaces
      for iz=1:Param.Grid.nz
        if y==Param.Grid.Faces(iF).Mid(2)
          cLoc=reshape(c(CG.Faces(iF).Glob,iz),OrdPoly+1,OrdPoly+1);
          dd=2/OrdPoly;
          eta0=-1;
          ksi0=-1;
          for iRef=1:OrdPoly
            ksi1=ksi0+dd;
            cc=0;
            for j=1:OrdPoly+1
              for i=1:OrdPoly+1
                cc=cc+Lagrange(0.5*(ksi0+ksi1),CG.xw,i)*Lagrange((eta0),CG.xw,j)...
                  *cLoc(i,j);
              end
            end
            [X(1,:)]=Trans(ksi0,eta0...
              ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
            [X(2,:)]=Trans(ksi1,eta0...
              ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
            [X(3,:)]=Trans(ksi1,eta0...
              ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
            [X(4,:)]=Trans(ksi0,eta0...
              ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
            fill3(X(:,1),X(:,2),X(:,3),cc);
            ksi0=ksi1;
          end
        end
      end
    end
  case 'YZ'
    x=Slice.x;
    for iF=1:Param.Grid.NumFaces
      for iz=1:Param.Grid.nz
        if x==Param.Grid.Faces(iF).Mid(1)
          cLoc=reshape(c(CG.Faces(iF).Glob,iz),OrdPoly+1,OrdPoly+1);
          dd=2/OrdPoly;
          eta0=-1;
          ksi0=-1;
          for jRef=1:OrdPoly
            eta1=eta0+dd;
            cc=0;
            for j=1:OrdPoly+1
              for i=1:OrdPoly+1
                cc=cc+Lagrange(ksi0,CG.xw,i)*Lagrange(0.5*(eta0+eta1),CG.xw,j)...
                  *cLoc(i,j);
              end
            end
            [X(1,:)]=Trans(ksi0,eta0...
              ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
            [X(2,:)]=Trans(ksi0,eta1...
              ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
            [X(3,:)]=Trans(ksi0,eta1...
              ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
            [X(4,:)]=Trans(ksi0,eta0...
              ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
            fill3(X(:,1),X(:,2),X(:,3),cc);
            eta0=eta1;
          end
        end
      end
    end
end
fig=fig+1;
hold off
end

