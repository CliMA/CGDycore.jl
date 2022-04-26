function vtkCG(c,CG,Global,vtkGrid,vtk)
OrdPoly=CG.OrdPoly;
OrdPrint=Global.Output.OrdPrint
nz=Global.Grid.nz;
NF=Global.Grid.NumFaces;
NumV=size(c,3);
ivtkc=0;
vtkc=zeros(nz*NF*OrdPrint*OrdPrint,NumV);
for iF=1:Global.Grid.NumFaces
  for iz=1:Global.Grid.nz
    dd=2/OrdPrint;
    eta0=-1;
    for jRef=1:OrdPrint
      ksi0=-1;
      eta1=eta0+dd;
      for iRef=1:OrdPrint
        ksi1=ksi0+dd;
        cc=zeros(NumV,1);
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            cc[:]=cc[:]+reshape(Lagrange(0.5*(ksi0+ksi1),CG.xw,i)*
              Lagrange(0.5*(eta0+eta1),CG.xw,j)*
              c[iz,CG.Glob[i,j,iF],:],NumV,1)
          end
        end
        ivtkc=ivtkc+1;
        vtkc[ivtkc,:]=cc[:];
        ksi0=ksi1;
      end
      eta0=eta1;
    end
  end
end

vtkS="$vtk"
vtkWriteHex(Global.Output.vtkFileName * vtkS * ".vtk",
  vtkGrid.vtkP,
  vtkGrid.ConnectivityList,vtkc,Global.Output.cNames)
vtk=vtk+1;
return vtk
end
