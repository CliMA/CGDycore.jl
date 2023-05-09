function InitSphere(OrdPoly,OrdPolyZ,nz,nPanel,H,GridType,TopoS,Decomp,stretch,Model,Phys)    

  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1 
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = InitParallelCom()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  Topography=(TopoS=TopoS,H=H,Rad=Phys.RadEarth)
  TimeStepper = InitTimeStepper()
  Grid = InitGrid(nz,Topography)

  if GridType == "HealPix"
  # Grid=CGDycore.InputGridH("Grid/mesh_H12_no_pp.nc",
  # CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
    Grid=InputGridH("Grid/mesh_H24_no_pp.nc", OrientFaceSphere,Phys.RadEarth,Grid)
  elseif GridType == "SQuadGen"
    Grid = InputGrid("Grid/baroclinic_wave_2deg_x4.g",OrientFaceSphere,Phys.RadEarth,Grid)
  elseif GridType == "CubedSphere"
    Grid = CubedGrid(nPanel,OrientFaceSphere,Phys.RadEarth,Grid)
  end

  if Decomp == "Hilbert"
    HilbertFaceSphere!(Grid)
    CellToProc = Decompose(Grid,ProcNumber)
  elseif Decomp == "EqualArea"
    CellToProc = DecomposeEqualArea(Grid,ProcNumber)
  else
    CellToProc = ones(Int,Grid.NumFaces)
    println(" False Decomp method ")
  end

  SubGrid = ConstructSubGrid(Grid,CellToProc,Proc)

  if stretch
    sigma = 1.0
    lambda = 3.16
#   AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
    AddExpStretchVerticalGrid!(SubGrid,nz,H,30.0,1500.0)
  else
    AddVerticalGrid!(SubGrid,nz,H)
  end

  Exchange = InitExchangeCG(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Model.HorLimit)
  Output = OutputStruct(Topography)
  Global = GlobalStruct{Float64}(SubGrid,Model,TimeStepper,ParallelCom,Phys,Output,Exchange,OrdPoly+1,nz,
    Model.NumV,Model.NumTr,())
  Global.Metric = Metric(OrdPoly+1,OrdPolyZ+1,SubGrid.NumFaces,nz)
  (CG,Global) = DiscretizationCG(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global)
  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = vtkInit3D(1,CGDycore.TransSphereX,CG,Global)
  unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp

  if TopoS == "EarthOrography"
    zS = CGDycore.Orography(CG,Global)
    Output.RadPrint = H
    Output.Flat=false
    nzTemp = Global.Grid.nz
    Global.Grid.nz = 1
    vtkCacheOrography = vtkInit2D(CG.OrdPoly,TransSphereX,CG,Global)
    unstructured_vtkOrography(zS,vtkCacheOrography, Global.Grid.NumFaces, CG,  Proc, ProcNumber)
    Global.Grid.nz = nzTemp
  end
  return Global
end  
