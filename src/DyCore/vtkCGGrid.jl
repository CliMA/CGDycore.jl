function [vtkGrid] = vtkCGGrid(CG,Trans,Topo,Param)
OrdPoly=CG.OrdPoly;
nz=Param.Grid.nz;
NF=Param.Grid.NumFaces;
if strcmp(Param.Grid.Form,'Sphere') && Param.Flat
  dTol=2*pi/max(Param.nPanel-1,1);
end
lam=zeros(8,1);
theta=zeros(8,1);

ivtkP=0;
vtkGrid.vtkP=zeros(3,8*Param.Grid.NumFaces*OrdPoly*OrdPoly*nz);
vtkGrid.ConnectivityList=reshape(1:1:8*NF*OrdPoly*OrdPoly*nz...
  ,8,NF*OrdPoly*OrdPoly*nz);
ivtkc=0;
for iF=1:Param.Grid.NumFaces
  for iz=1:Param.Grid.nz
    dd=2/OrdPoly;
    eta0=-1;
    for jRef=1:OrdPoly
      ksi0=-1;
      eta1=eta0+dd;
      for iRef=1:OrdPoly
        ksi1=ksi0+dd;
        [X(1,:)]=Trans(ksi0,eta0...
          ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
        [X(2,:)]=Trans(ksi1,eta0...
          ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
        [X(3,:)]=Trans(ksi1,eta1...
          ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
        [X(4,:)]=Trans(ksi0,eta1...
          ,Param.Grid.z(iz),Param.Grid.Faces(iF),Topo,Param);
        [X(5,:)]=Trans(ksi0,eta0...
          ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
        [X(6,:)]=Trans(ksi1,eta0...
          ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
        [X(7,:)]=Trans(ksi1,eta1...
          ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
        [X(8,:)]=Trans(ksi0,eta1...
          ,Param.Grid.z(iz+1),Param.Grid.Faces(iF),Topo,Param);
        if strcmp(Param.Grid.Form,'Sphere') && Param.Flat
          for i=1:8
          [lam(i),theta(i)]=cart2sphere(X(i,1),X(i,2),X(i,3));
          end
          
          lammin=min(lam);
          lammax=max(lam);
          if abs(lammin-lammax)>2*pi-dTol
            for i=1:8
              if lam(i)<pi
                lam(i)=lam(i)+2*pi;
                if lam(i)>3*pi
                  lam(i)=lam(i)-2*pi;
                end
              end
            end
          end
          for i=1:4
            ivtkP=ivtkP+1;
            vtkGrid.vtkP(:,ivtkP)=[lam(i),theta(i),Param.Grid.z(iz)/Param.H*3];
          end
          for i=5:8
            ivtkP=ivtkP+1;
            vtkGrid.vtkP(:,ivtkP)=[lam(i),theta(i),Param.Grid.z(iz+1)/Param.H*3];
          end
        else
          for i=1:8
            ivtkP=ivtkP+1;
            vtkGrid.vtkP(:,ivtkP)=[X(i,1),X(i,2),X(i,3)];
          end
        end
        ivtkc=ivtkc+1;
        ksi0=ksi1;
      end
      eta0=eta1;
    end
  end
end
end

