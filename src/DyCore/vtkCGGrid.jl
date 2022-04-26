mutable struct vtkCGGridStruct
  vtkP::Array{Float32, 2}
  ConnectivityList::Array{Int64, 2}
  OrdPrint::Int
end  

function vtkCGGrid(CG,Trans,Topo,Global)
OrdPoly=CG.OrdPoly;
OrdPrint=Global.Output.OrdPrint
nz=Global.Grid.nz;
NF=Global.Grid.NumFaces;
if Global.Grid.Form == "Sphere" && Global.Output.Flat
  dTol=2*pi/max(Global.Output.nPanel-1,1);
end
lam=zeros(8,1);
theta=zeros(8,1);
z=zeros(8,1);

ivtkP=0;
vtkP=zeros(Float32,3,8*NF*OrdPoly*OrdPoly*nz);
ConnectivityList=reshape(1:1:8*NF*OrdPoly*OrdPoly*nz,
  8,NF*OrdPoly*OrdPoly*nz);
ivtkc=0;
X = zeros(8,3)
for iF=1:Global.Grid.NumFaces
  for iz=1:Global.Grid.nz
    dd=2/OrdPrint;
    eta0=-1;
    for jRef=1:OrdPrint
      ksi0=-1;
      eta1=eta0+dd;
      for iRef=1:OrdPrint
        ksi1=ksi0+dd;
#       X[1,:]=Trans(ksi0,eta0, Global.Grid.z[iz],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[2,:]=Trans(ksi1,eta0, Global.Grid.z[iz],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[3,:]=Trans(ksi1,eta1, Global.Grid.z[iz],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[4,:]=Trans(ksi0,eta1, Global.Grid.z[iz],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[5,:]=Trans(ksi0,eta0, Global.Grid.z[iz+1],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[6,:]=Trans(ksi1,eta0, Global.Grid.z[iz+1],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[7,:]=Trans(ksi1,eta1, Global.Grid.z[iz+1],Global.Grid.Faces[iF],Topo,Global.Output);
#       X[8,:]=Trans(ksi0,eta1, Global.Grid.z[iz+1],Global.Grid.Faces[iF],Topo,Global.Output);
        X[1,:]=Trans(ksi0,eta0, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[2,:]=Trans(ksi1,eta0, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[3,:]=Trans(ksi1,eta1, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[4,:]=Trans(ksi0,eta1, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[5,:]=Trans(ksi0,eta0, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[6,:]=Trans(ksi1,eta0, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[7,:]=Trans(ksi1,eta1, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        X[8,:]=Trans(ksi0,eta1, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global);
        if Global.Grid.Form == "Sphere" && Global.Output.Flat
          for i=1:8
          (lam[i],theta[i],z[i])=cart2sphere(X[i,1],X[i,2],X[i,3]);
          end

          lammin=minimum(lam);
          lammax=maximum(lam);
          if abs(lammin-lammax)>2*pi-dTol
            for i=1:8
              if lam[i]<pi
                lam[i]=lam[i]+2*pi;
                if lam[i]>3*pi
                  lam[i]=lam[i]-2*pi;
                end
              end
            end
          end
          for i=1:8
            ivtkP=ivtkP+1;
            vtkP[:,ivtkP]=[lam[i],theta[i],max(z[i]-Global.Output.RadPrint,0.0)/Global.Output.H*3];
          end
#         for i=5:8
#           ivtkP=ivtkP+1;
#           vtkP[:,ivtkP]=[lam[i],theta[i],Global.Grid.z[iz+1]/Global.Output.H*3];
#         end
        else
          for i=1:8
            ivtkP=ivtkP+1;
            vtkP[:,ivtkP]=[X[i,1],X[i,2],X[i,3]];
          end
        end
        ivtkc=ivtkc+1;
        ksi0=ksi1;
      end
      eta0=eta1;
    end
  end
end  
return vtkCGGridStruct(
  vtkP,
  ConnectivityList,
  OrdPrint
  )
end

