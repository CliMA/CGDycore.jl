#function testNHHeldSuarezSphere
using SphericalGeometry
using CGDycore


OrdPoly = 4
OrdPolyZ = 1
nz = 2
nPanel = 4
NF = 6 * nPanel * nPanel

# Cache
cache=CGDycore.Cache(OrdPoly, OrdPolyZ, nz, NF)

# Physical parameters
Param=CGDycore.PhysParameters(cache);

# Grid
Param.nPanel=nPanel;
Param.H=30000;
Param.RadEarth=100
Param.Grid=CGDycore.CubedGrid(Param.nPanel,CGDycore.OrientFaceSphere,Param);

Param.Grid.nz=nz;
Param.Grid.zP=zeros(nz,1);
Param.Grid.z=zeros(nz+1,1);
Param.Grid.dz=Param.H/nz;
Param.Grid.zP[1]=Param.Grid.dz/2;
for i=2:nz
  Param.Grid.zP[i]=Param.Grid.zP[i-1]+Param.Grid.dz;
end
for i=2:nz+1
  Param.Grid.z[i]=Param.Grid.z[i-1]+Param.Grid.dz;
end


# Output
Param.RadPrint=Param.H;
Param.Flat=true;
#SphericalGrid
NumLon=20
NumLat=10
dlon=360.0/NumLon
dlat=180.0/NumLat
lon=zeros(NumLon,1)
lat=zeros(NumLat,1)
ksi=zeros(2,NumLon,NumLat)
Pos=zeros(Int,NumLon,NumLat)
lon1=pi/3
lat1=pi/10
CGDycore.FindPointInCell(lon1,lat1,Param.Grid)
lon[1]=0.5*dlon
for i=2:NumLon
  lon[i]=lon[i-1]+dlon  
end  
lat[1]=-pi/2+dlat/2
for j=2:NumLat
  lat[j]=lat[j-1]+dlat  
end  
for i=1:NumLat
  for j=1:NumLat
    (Pos[i,j],ksi[:,i,j])=CGDycore.FindPointInCell(lon[i],lat[j],Param.Grid)
  end
end  
println("iCell ",AAA)


# Initial conditions
U=zeros(CG.NumG,nz,Param.NumV);
U[:,:,Param.RhoPos]=CGDycore.Project(CGDycore.fRho,CG,Param);
(U[:,:,Param.uPos],U[:,:,Param.vPos])=CGDycore.ProjectVec(CGDycore.fVel,CG,Param);
U[:,:,Param.ThPos]=CGDycore.Project(CGDycore.fTheta,CG,Param).*U[:,:,Param.RhoPos];


# Integration
CFL=0.125;
dtau=600;
time=[0];

IntMethod="Rosenbrock";
# IntMethod="RungeKutta";
if strcmp(IntMethod,"Rosenbrock")
  dtau=600;
else
  dtau=8;
end
Param.RK=CGDycore.RungeKuttaMethod("RK4");
Param.ROS=CGDycore.RosenbrockMethod("SSP-Knoth");
SimDays=1000;
# SimDays=1;
PrintDay=10;
nIter=24*3600*SimDays/dtau;
PrintInt=24*3600*PrintDay/dtau;
# Print initial conditions
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
#
str = IntMethod
OP=CG.OrdPoly+1;
NF=Param.Grid.NumFaces;
nz=Param.Grid.nz;
Param.CacheF1=zeros(OP,OP,NF,nz+1);
Param.CacheF2=zeros(OP,OP,NF,nz+1);
Param.CacheF3=zeros(OP,OP,NF,nz+1);
Param.CacheF4=zeros(OP,OP,NF,nz+1);
Param.CacheF5=zeros(OP,OP,NF,nz+1);
Param.CacheF6=zeros(OP,OP,NF,nz+1);
Param.CacheC1 = view(Param.CacheF1,:,:,:,1:nz)
Param.CacheC2 = view(Param.CacheF2,:,:,:,1:nz)
Param.CacheC3 = view(Param.CacheF3,:,:,:,1:nz)
Param.CacheC4 = view(Param.CacheF4,:,:,:,1:nz)
Param.CacheC5 = view(Param.CacheF5,:,:,:,1:nz)
Param.CacheC6 = view(Param.CacheF6,:,:,:,1:nz)
Param.Cache1=zeros(CG.NumG,nz)
Param.Cache2=zeros(CG.NumG,nz)
Param.Cache3=zeros(CG.NumG,nz)
Param.Cache4=zeros(CG.NumG,nz)
Param.Pres=zeros(OP,OP,NF,nz)
Param.KE=zeros(OP,OP,NF,nz)
Param.FCG=zeros(OP,OP,NF,nz,size(U,3))
Param.k=zeros(size(U)..., Param.ROS.nStage);
Param.fV=zeros(size(U))
Param.Vn=zeros(size(U))
Param.RhoCG=zeros(OP,OP,NF,nz)
Param.v1CG=zeros(OP,OP,NF,nz)
Param.v2CG=zeros(OP,OP,NF,nz)
Param.wCG=zeros(OP,OP,NF,nz+1)
Param.wCCG=zeros(OP,OP,NF,nz+1)
Param.ThCG=zeros(OP,OP,NF,nz)
Param.J = CGDycore.JacStruct(CG.NumG,nz)
if str == "Rosenbrock"
    @time begin
      for i=1:nIter
        Δt = @elapsed begin
          CGDycore.RosenbrockSchur!(U,dtau,CGDycore.FcnNHCurlVec!,CGDycore.JacSchur,CG,Param);
          time[1] += dtau;
          if mod(i,PrintInt)==0
            Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
          end
        end
        percent = i/nIter*100
        @info "Iteration: $i took $Δt, $percent% complete"
      end
    end

elseif str == "RungeKutta"
    for i=1:nIter
      @info "Iteration: $i"

      U .= CGDycore.RungeKuttaExplicit(U,dtau,CGDycore.FcnNHCurlVec,CG,Param);

      time[1] += dtau;
      if mod(i,PrintInt)==0
        Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param);
      end
    end
else
  error("Bad str")
end
Param.vtk=CGDycore.vtkOutput(U,vtkGrid,CG,Param)

@info "Success!"
