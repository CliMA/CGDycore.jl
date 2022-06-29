using CGDycore

function testDecompose()

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



ProcNumber = 3

Proc = 1
Grid1=CGDycore.Grid(nz,Topography)
Grid1=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid1)
CellToProc = CGDycore.Decompose(Grid1,ProcNumber)
SubGrid1 = CGDycore.ConstructSubGrid(Grid1,CellToProc,Proc)
CGDycore.AddVerticalGrid!(SubGrid1,nz,H)
Output1=CGDycore.Output(Topography)
Global1 = CGDycore.Global(SubGrid1,Model,Phys,Output1,OrdPoly+1,nz,NumV,NumTr,())
Global1.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid1.NumFaces,nz)
(CG1,Global1)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global1)
# Output
  Output1.vtkFileName="BaroWaveMoist1"
  Output1.vtk=0
  Output1.Flat=false
  Output1.nPanel=nPanel
  Output1.RadPrint=H
  Output1.H=H
  Output1.cNames = [
    "u",
]

Proc = 2
Grid2=CGDycore.Grid(nz,Topography)
Grid2=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid2)
CellToProc = CGDycore.Decompose(Grid2,ProcNumber)
SubGrid2 = CGDycore.ConstructSubGrid(Grid2,CellToProc,Proc)
CGDycore.AddVerticalGrid!(SubGrid2,nz,H)
Output2=CGDycore.Output(Topography)
Global2 = CGDycore.Global(SubGrid2,Model,Phys,Output2,OrdPoly+1,nz,NumV,NumTr,())
Global2.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid2.NumFaces,nz)
(CG2,Global2)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global2)
# Output
  Output2.vtkFileName="BaroWaveMoist2"
  Output2.vtk=0
  Output2.Flat=false
  Output2.nPanel=nPanel
  Output2.RadPrint=H
  Output2.H=H
  Output2.cNames = [
    "u",
]

Proc = 3
Grid3=CGDycore.Grid(nz,Topography)
Grid3=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid3)
CellToProc = CGDycore.Decompose(Grid3,ProcNumber)
SubGrid3 = CGDycore.ConstructSubGrid(Grid3,CellToProc,Proc)
CGDycore.AddVerticalGrid!(SubGrid3,nz,H)
Output3=CGDycore.Output(Topography)
Global3 = CGDycore.Global(SubGrid3,Model,Phys,Output3,OrdPoly+1,nz,NumV,NumTr,())
Global3.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid3.NumFaces,nz)
(CG3,Global3)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global3)
# Output
  Output3.vtkFileName="BaroWaveMoist3"
  Output3.vtk=0
  Output3.Flat=false
  Output3.nPanel=nPanel
  Output3.RadPrint=H
  Output3.H=H
  Output3.cNames = [
    "u",
]

# Output
  Output1.OrdPrint=CG1.OrdPoly
  vtkGrid1=CGDycore.vtkCGGrid(CG1,CGDycore.TransSphereX,CGDycore.Topo,Global1)

  Output2.OrdPrint=CG2.OrdPoly
  vtkGrid2=CGDycore.vtkCGGrid(CG2,CGDycore.TransSphereX,CGDycore.Topo,Global2)

  Output3.OrdPrint=CG3.OrdPoly
  vtkGrid3=CGDycore.vtkCGGrid(CG3,CGDycore.TransSphereX,CGDycore.Topo,Global3)

  U1=zeros(nz,CG1.NumG,Model.NumV+Model.NumTr)
  U1[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG1,Global1)
  (U1[:,:,Model.uPos],U1[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG1,Global1)
  U1[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG1,Global1).*U1[:,:,Model.RhoPos]
  if NumTr>0
    U1[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG1,Global1).*U1[:,:,Model.RhoPos]
  end   
  U1[:,:,2] .= 1.0

  U2=zeros(nz,CG2.NumG,Model.NumV+Model.NumTr)
  U2[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG2,Global2)
  (U2[:,:,Model.uPos],U2[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG2,Global2)
  U2[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG2,Global2).*U2[:,:,Model.RhoPos]
  if NumTr>0
    U2[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG2,Global2).*U2[:,:,Model.RhoPos]
  end
  U2[:,:,2] .= 2.0

  U3=zeros(nz,CG3.NumG,Model.NumV+Model.NumTr)
  U3[:,:,Model.RhoPos]=CGDycore.Project(CGDycore.fRho,0.0,CG3,Global3)
  (U3[:,:,Model.uPos],U3[:,:,Model.vPos])=CGDycore.ProjectVec(CGDycore.fVel,0.0,CG3,Global3)
  U3[:,:,Model.ThPos]=CGDycore.Project(CGDycore.fTheta,0.0,CG3,Global3).*U3[:,:,Model.RhoPos]
  if NumTr>0
    U3[:,:,Model.RhoVPos+Model.NumV]=CGDycore.Project(CGDycore.fQv,0.0,CG3,Global3).*U3[:,:,Model.RhoPos]
  end
  U3[:,:,2] .= 3.0

  DictE1=Dict()
  for iE = 1 : SubGrid1.NumInBoundEdges
    DictE1[SubGrid1.Edges[SubGrid1.InBoundEdges[iE]].EG] = (SubGrid1.Edges[SubGrid1.InBoundEdges[iE]].E,
      SubGrid1.InBoundEdgesP[iE])
  end  
  DictE2=Dict()
  for iE = 1 : SubGrid2.NumInBoundEdges
    DictE2[SubGrid2.Edges[SubGrid2.InBoundEdges[iE]].EG] = (SubGrid2.Edges[SubGrid2.InBoundEdges[iE]].E,
      SubGrid2.InBoundEdgesP[iE])
  end  
  DictE3=Dict()
  for iE = 1 : SubGrid3.NumInBoundEdges
    DictE3[SubGrid3.Edges[SubGrid3.InBoundEdges[iE]].EG] = (SubGrid3.Edges[SubGrid3.InBoundEdges[iE]].E,
      SubGrid3.InBoundEdgesP[iE])
  end  

  #Proc1
   NumNeiProc1 = 2 
   NeiProc1 = [2, 3]
   GlobBuffer1 = Dict()
   LocBuffer1 = Dict()
   for iP = 1 : NumNeiProc1
     LocTemp=zeros(Int,0)  
     GlobTemp=zeros(Int,0)  
     for iEB = 1 : SubGrid1.NumInBoundEdges
       if SubGrid1.InBoundEdgesP[iEB] == NeiProc1[iP]
         iE = SubGrid1.InBoundEdges[iEB] 
         push!(GlobTemp,SubGrid1.Edges[iE].EG)
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid1.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end  
     end
     LocBuffer1[NeiProc1[iP]] = LocTemp
     GlobBuffer1[NeiProc1[iP]] = GlobTemp
   end  

  #Proc2
   NumNeiProc2 = 2 
   NeiProc2 = [1, 3]
   LocBuffer2 = Dict()
   GlobBuffer2 = Dict()
   for iP = 1 : NumNeiProc2
     LocTemp=zeros(Int,0)  
     GlobTemp=zeros(Int,0)  
     for iEB = 1 : SubGrid2.NumInBoundEdges
       if SubGrid2.InBoundEdgesP[iEB] == NeiProc2[iP]
         iE = SubGrid2.InBoundEdges[iEB] 
         push!(GlobTemp,SubGrid2.Edges[iE].EG)
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid2.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end  
     end
     LocBuffer2[NeiProc2[iP]] = LocTemp
     GlobBuffer2[NeiProc2[iP]] = GlobTemp
   end  

  #Proc3
   NumNeiProc3 = 2 
   NeiProc3 = [1, 2]
   LocBuffer3 = Dict()
   GlobBuffer3 = Dict()
   for iP = 1 : NumNeiProc3
     LocTemp=zeros(Int,0)  
     GlobTemp=zeros(Int,0)  
     for iEB = 1 : SubGrid3.NumInBoundEdges
       if SubGrid3.InBoundEdgesP[iEB] == NeiProc3[iP]
         iE = SubGrid3.InBoundEdges[iEB] 
         push!(GlobTemp,SubGrid3.Edges[iE].EG)
         for k = 1 : OrdPoly - 1
           push!(LocTemp,k + SubGrid3.NumNodes + (iE - 1)*(OrdPoly - 1))
         end
       end  
     end
     LocBuffer3[NeiProc3[iP]] = LocTemp
     GlobBuffer3[NeiProc3[iP]] = GlobTemp
   end  


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
end
