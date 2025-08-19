import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore, FEMSei, FiniteVolumes
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using LinearAlgebra

# Model
parsed_args = DyCore.parse_commandline()
Problem = parsed_args["Problem"]
ProfRho = parsed_args["ProfRho"]
ProfTheta = parsed_args["ProfTheta"]
PertTh = parsed_args["PertTh"]
ProfVel = parsed_args["ProfVel"]
ProfVelGeo = parsed_args["ProfVelGeo"]
RhoVPos = parsed_args["RhoVPos"]
RhoCPos = parsed_args["RhoCPos"]
RhoIPos = parsed_args["RhoIPos"]
RhoRPos = parsed_args["RhoRPos"]
HorLimit = parsed_args["HorLimit"]
Upwind = parsed_args["Upwind"]
Damping = parsed_args["Damping"]
Relax = parsed_args["Relax"]
StrideDamp = parsed_args["StrideDamp"]
Geos = parsed_args["Geos"]
Coriolis = parsed_args["Coriolis"]
CoriolisType = parsed_args["CoriolisType"]
Buoyancy = parsed_args["Buoyancy"]
Equation = parsed_args["Equation"]
RefProfile = parsed_args["RefProfile"]
ProfpBGrd = parsed_args["ProfpBGrd"]
ProfRhoBGrd = parsed_args["ProfRhoBGrd"]
Microphysics = parsed_args["Microphysics"]
TypeMicrophysics = parsed_args["TypeMicrophysics"]
RelCloud = parsed_args["RelCloud"]
Rain = parsed_args["Rain"]
Source = parsed_args["Source"]
Forcing = parsed_args["Forcing"]
VerticalDiffusion = parsed_args["VerticalDiffusion"]
JacVerticalDiffusion = parsed_args["JacVerticalDiffusion"]
JacVerticalAdvection = parsed_args["JacVerticalAdvection"]
SurfaceFlux = parsed_args["SurfaceFlux"]
SurfaceFluxMom = parsed_args["SurfaceFluxMom"]
NumV = parsed_args["NumV"]
NumTr = parsed_args["NumTr"]
Curl = parsed_args["Curl"]
ModelType = parsed_args["ModelType"]
Thermo = parsed_args["Thermo"]
# Parallel
Decomp = parsed_args["Decomp"]
# Time integration
SimDays = parsed_args["SimDays"]
SimHours = parsed_args["SimHours"]
SimMinutes = parsed_args["SimMinutes"]
SimSeconds = parsed_args["SimSeconds"]
StartAverageDays = parsed_args["StartAverageDays"]
dtau = parsed_args["dtau"]
IntMethod = parsed_args["IntMethod"]
Table = parsed_args["Table"]
# Grid
nz = parsed_args["nz"]
nPanel = parsed_args["nPanel"]
H = parsed_args["H"]
Stretch = parsed_args["Stretch"]
StretchType = parsed_args["StretchType"]
TopoS = parsed_args["TopoS"]
GridType = parsed_args["GridType"]
RadEarth = parsed_args["RadEarth"]
# CG Element
OrdPoly = parsed_args["OrdPoly"]
# Viscosity
HyperVisc = parsed_args["HyperVisc"]
HyperDCurl = parsed_args["HyperDCurl"]
HyperDGrad = parsed_args["HyperDGrad"]
HyperDRhoDiv = parsed_args["HyperDRhoDiv"]
HyperDDiv = parsed_args["HyperDDiv"]
HyperDDivW = parsed_args["HyperDDivW"]
# Output
PrintDays = parsed_args["PrintDays"]
PrintHours = parsed_args["PrintHours"]
PrintMinutes = parsed_args["PrintMinutes"]
PrintSeconds = parsed_args["PrintSeconds"]
PrintStartTime = parsed_args["PrintStartTime"]
Flat = parsed_args["Flat"]

# Device
Device = parsed_args["Device"]
GPUType = parsed_args["GPUType"]
FloatTypeBackend = parsed_args["FloatTypeBackend"]
NumberThreadGPU = parsed_args["NumberThreadGPU"]

MPI.Init()

Device = "CPU"
FloatTypeBackend = "Float64"

if Device == "CPU" 
  backend = CPU()
elseif Device == "GPU" 
  if GPUType == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(true)
#   CUDA.device!(MPI.Comm_rank(MPI.COMM_WORLD))
  elseif GPUType == "AMD"
    backend = ROCBackend()
    AMDGPU.allowscalar(false)
  elseif GPUType == "Metal"
    backend = MetalBackend()
    Metal.allowscalar(true)
  end
else
  backend = CPU()
end

if FloatTypeBackend == "Float64"
  FTB = Float64
elseif FloatTypeBackend == "Float32"
  FTB = Float32
else
  @show "False FloatTypeBackend"
  stop
end

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

#ModelParameters
Model = DyCore.ModelStruct{FTB}()

RefineLevel = 5
RadEarth = 1.0
nz = 1
nPanel = 40
OP = 4
DoF = OP * OP
OPZ = 2
Decomp = ""
Decomp = "EqualArea"

GridType = "CubedSphere"
Grid, Exchange = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,
  Model,ParallelCom;order=false)
  iF = 1000
  iF = 333
  F = zeros(4,3)
  F[1,1] = Grid.Faces[iF].P[1].x
  F[1,2] = Grid.Faces[iF].P[1].y
  F[1,3] = Grid.Faces[iF].P[1].z
  F[2,1] = Grid.Faces[iF].P[2].x
  F[2,2] = Grid.Faces[iF].P[2].y
  F[2,3] = Grid.Faces[iF].P[2].z
  F[3,1] = Grid.Faces[iF].P[3].x
  F[3,2] = Grid.Faces[iF].P[3].y
  F[3,3] = Grid.Faces[iF].P[3].z
  F[4,1] = Grid.Faces[iF].P[4].x
  F[4,2] = Grid.Faces[iF].P[4].y
  F[4,3] = Grid.Faces[iF].P[4].z


  Rad = 1.e5
  zero = eltype(F)(0)
  one = eltype(F)(1)
  half = eltype(F)(1/2)
  quarter = eltype(F)(1/4)
  ksi1 = 0.75
  ksi2 = 0.75
  dXdx = zeros(2,2)
  X1 = quarter * (F[1,1] * (one-ksi1)*(one-ksi2) +
   F[2,1] * (one+ksi1)*(one-ksi2) +
   F[3,1] * (one+ksi1)*(one+ksi2) +
   F[4,1] * (one-ksi1)*(one+ksi2))
  X2 = quarter * (F[1,2] * (one-ksi1)*(one-ksi2) +
   F[2,2] * (one+ksi1)*(one-ksi2) +
   F[3,2] * (one+ksi1)*(one+ksi2) +
   F[4,2] * (one-ksi1)*(one+ksi2))
  X3 = quarter * (F[1,3] * (one-ksi1)*(one-ksi2) +
   F[2,3] * (one+ksi1)*(one-ksi2) +
   F[3,3] * (one+ksi1)*(one+ksi2) +
   F[4,3] * (one-ksi1)*(one+ksi2))

  r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)
  f = Rad / r
  X1 = X1 / r
  X2 = X2 / r
  X3 = X3 / r
  (lam,theta)=Grids.cart2sphere(X1,X2,X3)

  DD=@SArray([-sin(lam) cos(lam) zero;
      zero       zero     one])

  sinlam = sin(lam)
  coslam = cos(lam)
  sinth = sin(theta)
  costh = cos(theta)
  a11 = sinlam * sinlam * costh * costh + sinth * sinth
  a12 = -sinlam * coslam * costh * costh
  a13 = -coslam * sinth * costh
  a21 = a12
  a22 = coslam * coslam * costh * costh + sinth * sinth
  a23 = -sinlam * sinth * costh
  a31 = -coslam * sinth
  a32 = -sinlam * sinth
  a33 = costh
  A = @SArray([a11 a12 a13;
      a21 a22 a23;
      a31 a32 a33])

  B = @SArray([F[1,1] F[2,1] F[3,1] F[4,1];
       F[1,2] F[2,2] F[3,2] F[4,2];
       F[1,3] F[2,3] F[3,3] F[4,3]])

  C = @SArray([-one+ksi2  -one+ksi1;
              one-ksi2  -one-ksi1;
              one+ksi2   one+ksi1;
             -one-ksi2   one-ksi1])
  D = quarter * f * DD * A * B * C
  dXdx[1,1] = D[1,1]	    
  dXdx[1,2] = D[1,2]	    
  dXdx[2,1] = D[2,1]	    
  dXdx[2,2] = D[2,2]	    

  detdXdx = det(dXdx) 
  dXdxI = inv(dXdx) * detdXdx


  P1 = Grid.Faces[iF].P[1]
  P2 = Grid.Faces[iF].P[2]
  P3 = Grid.Faces[iF].P[3]
  P4 = Grid.Faces[iF].P[4]

  XT1 =  0.25*(P1.x*(1-ksi1)*(1-ksi2)+
            P2.x*(1+ksi1)*(1-ksi2)+
            P3.x*(1+ksi1)*(1+ksi2)+
            P4.x*(1-ksi1)*(1+ksi2))
    
  XT2 =  0.25*(P1.y*(1-ksi1)*(1-ksi2)+
            P2.y*(1+ksi1)*(1-ksi2)+
            P3.y*(1+ksi1)*(1+ksi2)+
            P4.y*(1-ksi1)*(1+ksi2))
           
  XT3 =  0.25*(P1.z*(1-ksi1)*(1-ksi2)+
           P2.z*(1+ksi1)*(1-ksi2)+
            P3.z*(1+ksi1)*(1+ksi2)+
            P4.z*(1-ksi1)*(1+ksi2))
    
  X = SVector{3}(XT1,XT2,XT3)
  JP = @SArray[P1.x P2.x P3.x P4.x;
               P1.y P2.y P3.y P4.y;
               P1.z P2.z P3.z P4.z]

  J3 = @SArray([-0.25 + 0.25*ksi2  -0.25 + 0.25*ksi1
                 0.25 - 0.25*ksi2  -0.25 - 0.25*ksi1
                 0.25 + 0.25*ksi2   0.25 + 0.25*ksi1
                -0.25 - 0.25*ksi2   0.25 - 0.25*ksi1])
# X[i] = xT[i]/norm(XT)*(R + 0.5*(1-xi3)*h1 + (1+xi3)*h2)

  f1 = 1.0 / (XT1^2 + XT2^2 + XT3^2)^(1/2)
  f = Rad *(XT1^2 + XT2^2 + XT3^2)^(-3/2)
  dX1dXT1 = f * (XT2^2 + XT3^2)
  dX1dXT2= -f * XT1 * XT2 
  dX1dXT3= -f * XT1 * XT3
  dX2dXT1 = dX1dXT2 
  dX2dXT2 = f * (XT1^2+XT3^2)
  dX2dXT3 = -f * XT2 * XT3
  dX3dXT1 = dX1dXT3 
  dX3dXT2 = dX2dXT3 
  dX3dXT3 = f*(XT1^2+XT2^2)

  J1 = @SArray([dX1dXT1    dX2dXT1     dX3dXT1   
                dX1dXT2     dX2dXT2     dX3dXT2
                dX1dXT3     dX2dXT3     dX3dXT3])   
  J = J1 * JP * J3      

  detJ = norm(cross(J[:,1],J[:,2]))
  pinvJ = detJ * (J / (J' * J))

  (lon,lat)=Grids.cart2sphere(XT1,XT2,XT3)
  Rotate = zeros(3,3)
  Rotate[1,1] =           -sin(lon)
  Rotate[2,1] =-sin(lat) * cos(lon)
  Rotate[3,1] = cos(lat) * cos(lon)
  Rotate[1,2] =           cos(lon)
  Rotate[2,2] =-sin(lat) * sin(lon)
  Rotate[3,2] = cos(lat) * sin(lon)
  Rotate[1,3] =          0.0
  Rotate[2,3] = cos(lat)
  Rotate[3,3] = sin(lat)
  JJ = [J[:,1] J[:,2] [XT1 / f1 ;XT2 / f1 ;XT3 / f1 ]]

