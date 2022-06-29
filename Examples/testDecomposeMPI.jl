using CGDycore
using MPI

#function testDecompose()
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
print("$rank: \n")
size = MPI.Comm_size(comm)
ProcNumber = size
Proc = rank + 1
print("$Proc: \n")
print("$ProcNumber: \n")

OrdPoly = 4
nz = 1

OrdPolyZ=1
nPanel = 8
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 2


# Physical parameters
Phys=CGDycore.PhysParameters()

#ModelParameters
day=3600*24
Param=(T0E=310.0,
       T0P=240.0,
       B=2.0,
       K=3.0,
       LapseRate=0.005,
       U0=-0.5,
       PertR=1.0/6.0,
       Up=1.0,
       PertExpR=0.1,
       PertLon=pi/9.0,
       PertLat=2.0*pi/9.0,
       PertZ=15000.0,
       NBr=1.e-2,
       DeltaT=1,
       ExpDist=5,
       T0=300,
       T_init = 315,
       lapse_rate = -0.008,
       Deep=false,
       k_a=1/(40 * day),
       k_f=1/day,
       k_s=1/(4 * day),
       pert = 0.1,
       uMax = 1.0,
       vMax = 0.0,
       DeltaT_y=0,
       DeltaTh_z=-5,
       T_equator=315,
       T_min=200,
       sigma_b=7/10,
       z_D=20.0e3,
#      Boundary layer       
       C_E = 0.0044,
       p_pbl = 85000.0,
       p_strato = 10000.0,
#      Surface flux       
       CTr = 0.004,
#      Moist
       q_0 = 0.018,                # Maximum specific humidity (default: 0.018)
       q_t = 1e-12,
       )
Model = CGDycore.Model(Param)
# Initial conditions
  Model.Equation="CompressibleMoist"
  Model.NumV=NumV
  Model.NumTr=NumTr
  Model.Problem="BaroWaveMoistSphere"
  Model.ProfRho="BaroWaveSphere"
  Model.ProfTheta="BaroWaveMoistSphere"
  Model.ProfVel="BaroWaveSphere"
  Model.RhoPos=1
  Model.uPos=2
  Model.vPos=3
  Model.wPos=4
  Model.ThPos=5
  Model.RhoVPos = 1
  Model.RhoCPos = 2
  Model.HorLimit = false
  Model.Source = false
  Model.Upwind = true
  Model.Microphysics = false
  Model.RelCloud = 0.01
  Model.Damping = false
  Model.StrideDamp=20000.0
  Model.Relax = 1.0/100.0

# Grid
H = 30000.0
#H = 45000.0
Topography=(TopoS="EarthOrography",H=H,Rad=Phys.RadEarth)




Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
CellToProc = CGDycore.Decompose(Grid,ProcNumber)
SubGrid = CGDycore.ConstructSubGrid(Grid,CellToProc,Proc)
CGDycore.AddVerticalGrid!(SubGrid,nz,H)
Output=CGDycore.Output(Topography)
Global = CGDycore.Global(SubGrid,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())
Global.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
(CG,Global)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global)
# Output
  Output.vtkFileName=string("BaroWaveMoist",string(rank))
  Output.vtk=0
  Output.Flat=false
  Output.nPanel=nPanel
  Output.RadPrint=H
  Output.H=H
  Output.cNames = [
    "u",
]
@show Output.vtkFileName

# Output
  Output.OrdPrint=CG.OrdPoly
  vtkGrid=CGDycore.vtkCGGrid(CG,CGDycore.TransSphereX,CGDycore.Topo,Global)


  U = zeros(nz,CG.NumG,Model.NumV+Model.NumTr)
  U[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG,Global)
  (U[:,:,Model.uPos],U[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG,Global)
  U[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG,Global).*U[:,:,Model.RhoPos]
  if NumTr>0
    U[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG,Global).*U[:,:,Model.RhoPos]
  end   
  U[:,:,2] .= Proc

  DictE=Dict()
  for iE = 1 : SubGrid.NumInBoundEdges
    DictE[SubGrid.Edges[SubGrid.InBoundEdges[iE]].EG] = (SubGrid.Edges[SubGrid.InBoundEdges[iE]].E,
      SubGrid.InBoundEdgesP[iE])
  end  

   NumNeiProc = SubGrid.NumNeiProc
   NeiProc = SubGrid.NeiProc
   GlobBuffer = Dict()
   LocBuffer = Dict()
   for iP = 1 : NumNeiProc
     LocTemp=zeros(Int,0)  
     GlobTemp=zeros(Int,0)  
     for iEB = 1 : SubGrid.NumInBoundEdges
       if SubGrid.InBoundEdgesP[iEB] == NeiProc[iP]
         iE = SubGrid.InBoundEdges[iEB] 
         push!(GlobTemp,SubGrid.Edges[iE].EG)
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end  
     end
     LocBuffer[NeiProc[iP]] = LocTemp
     GlobBuffer[NeiProc[iP]] = GlobTemp
   end  

   GlobGetBuffer = Dict()
   rreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
   for iP = 1 : NumNeiProc
     GlobGetBuffer[NeiProc[iP]] = similar(GlobBuffer[NeiProc[iP]])
     tag = Proc + ProcNumber*NeiProc[iP]
     rreq[iP] = MPI.Irecv!(GlobGetBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, comm)
   end  
   sreq = MPI.Request[MPI.REQUEST_NULL for _ in (NeiProc .- 1)]
   for iP = 1 : NumNeiProc
     tag = NeiProc[iP] + ProcNumber*Proc
     sreq[iP] = MPI.Isend(GlobBuffer[NeiProc[iP]], NeiProc[iP] - 1, tag, comm)
   end  

   stats = MPI.Waitall!(rreq)
   stats = MPI.Waitall!(sreq)

   MPI.Barrier(comm)
#=

   GlobGetBuffer1 = Dict()
   GlobGetBuffer1[2] = GlobBuffer2[1]
   GlobGetBuffer1[3] = GlobBuffer3[1]
   GlobGetBuffer2 = Dict()
   GlobGetBuffer2[1] = GlobBuffer1[2]
   GlobGetBuffer2[3] = GlobBuffer3[2]
   GlobGetBuffer3 = Dict()
   GlobGetBuffer3[1] = GlobBuffer1[3]
   GlobGetBuffer3[2] = GlobBuffer2[3]


   #Proc1
   GetBuffer1 = Dict()
   for iP = 1 : NumNeiProc1
     GlobInd = GlobGetBuffer1[NeiProc1[iP]]  
     LocTemp=zeros(Int,0)
     iEB1 = 0
     for iEB = 1 : SubGrid1.NumInBoundEdges
       if SubGrid1.InBoundEdgesP[iEB] == NeiProc1[iP]
         iEB1 += 1  
         (iE,) = DictE1[GlobInd[iEB1]]  
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid1.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end
     end
     GetBuffer1[NeiProc1[iP]] = LocTemp
   end

   #Proc2
   GetBuffer2 = Dict()
   for iP = 1 : NumNeiProc2
     GlobInd = GlobGetBuffer2[NeiProc2[iP]]  
     LocTemp=zeros(Int,0)
     iEB2 = 0
     for iEB = 1 : SubGrid2.NumInBoundEdges
       if SubGrid2.InBoundEdgesP[iEB] == NeiProc2[iP]
         iEB2 += 1 
         (iE,) = DictE2[GlobInd[iEB2]]  
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid1.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end
     end
     GetBuffer2[NeiProc2[iP]] = LocTemp
   end

   #Proc3
   GetBuffer3 = Dict()
   for iP = 1 : NumNeiProc3
     GlobInd = GlobGetBuffer3[NeiProc3[iP]]  
     LocTemp=zeros(Int,0)
     iEB3 = 0
     for iEB = 1 : SubGrid3.NumInBoundEdges
       if SubGrid3.InBoundEdgesP[iEB] == NeiProc3[iP]
         iEB3 += 1   
         (iE,) = DictE3[GlobInd[iEB3]]  
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid1.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end
     end
     GetBuffer3[NeiProc3[iP]] = LocTemp
   end

  #Send
  #Proc1
   SendBuffer1 = Dict()
   for iP = 1 : NumNeiProc1
     SendBuffer1[NeiProc1[iP]] = U1[1,LocBuffer1[NeiProc1[iP]],2]
   end
  #Proc2
   SendBuffer2 = Dict()
   for iP = 1 : NumNeiProc2
     SendBuffer2[NeiProc2[iP]] = U2[1,LocBuffer2[NeiProc2[iP]],2]
   end
  #Proc3
   SendBuffer3 = Dict()
   for iP = 1 : NumNeiProc3
     SendBuffer3[NeiProc3[iP]] = U3[1,LocBuffer3[NeiProc3[iP]],2]
   end

  #Exchange 
   RecvBuffer1 = Dict()
   RecvBuffer1[2] = SendBuffer2[1]
   RecvBuffer1[3] = SendBuffer3[1]
   RecvBuffer2 = Dict()
   RecvBuffer2[1] = SendBuffer1[2]
   RecvBuffer2[3] = SendBuffer3[2]
   RecvBuffer3 = Dict()
   RecvBuffer3[1] = SendBuffer1[3]
   RecvBuffer3[2] = SendBuffer2[3]

  #Receive 
  #Proc1
   for iP = 1 : NumNeiProc1
     U1[1,GetBuffer1[NeiProc1[iP]],2] .+= RecvBuffer1[NeiProc1[iP]]
   end
  #Proc2
   for iP = 1 : NumNeiProc2
     U2[1,GetBuffer2[NeiProc2[iP]],2] .+= RecvBuffer2[NeiProc2[iP]]
   end
  #Proc3
   for iP = 1 : NumNeiProc3
     U3[1,GetBuffer3[NeiProc3[iP]],2] .+= RecvBuffer3[NeiProc3[iP]]
   end

  Global1.Output.vtk=CGDycore.vtkOutput(U1,vtkGrid1,CG1,Global1)
  Global2.Output.vtk=CGDycore.vtkOutput(U2,vtkGrid2,CG2,Global2)
  Global3.Output.vtk=CGDycore.vtkOutput(U3,vtkGrid3,CG3,Global3)
  =#
#end
