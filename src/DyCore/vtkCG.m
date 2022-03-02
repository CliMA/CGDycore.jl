function [vtk] = vtkCG(c,CG,Param,vtkGrid,vtk)
OrdPoly=CG.OrdPoly;
nz=Param.Grid.nz;
NF=Param.Grid.NumFaces;
NumV=size(c,3);
ivtkc=0;
vtkc=zeros(nz*NF*OrdPoly*OrdPoly,NumV);
for iF=1:Param.Grid.NumFaces
  for iz=1:Param.Grid.nz
    cLoc=reshape(c(CG.Faces(iF).Glob,iz,:),OrdPoly+1,OrdPoly+1,NumV);
    dd=2/OrdPoly;
    eta0=-1;
    for jRef=1:OrdPoly
      ksi0=-1;
      eta1=eta0+dd;
      for iRef=1:OrdPoly
        ksi1=ksi0+dd;
        cc=zeros(NumV,1);
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            cc(:)=cc(:)+reshape(Lagrange(0.5*(ksi0+ksi1),CG.xw,i)...
              *Lagrange(0.5*(eta0+eta1),CG.xw,j)...
              *cLoc(i,j,:),NumV,1);
          end
        end
        ivtkc=ivtkc+1;
        vtkc(ivtkc,:)=cc(:);
        ksi0=ksi1;
      end
      eta0=eta1;
    end
  end
end

vtkS=num2str(vtk);
vtkWriteHex(strcat(Param.vtkFileName,vtkS,'.vtk'),vtkGrid.vtkP...
  ,vtkGrid.ConnectivityList,vtkc,Param.cNames)
vtk=vtk+1;
end

