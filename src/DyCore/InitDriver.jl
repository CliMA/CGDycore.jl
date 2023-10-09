function InitSphere(backend,FT,OrdPoly,OrdPolyZ,nz,nPanel,H,GridType,Topography,Decomp,Model,Phys)    

  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1 
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = ParallelComStruct()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)
  Grid = GridStruct(nz,Topography)

  if GridType == "HealPix"
  # Grid=CGDycore.InputGridH("Grid/mesh_H12_no_pp.nc",
  # CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
    Grid=InputGridH("Grid/mesh_H24_no_pp.nc", OrientFaceSphere,Phys.RadEarth,Grid)
  elseif GridType == "SQuadGen"
    Grid = InputGrid("Grid/baroclinic_wave_2deg_x4.g",OrientFaceSphere,Phys.RadEarth,Grid)
  elseif GridType == "Msh"
    Grid = InputGridMsh("Grid/Quad.msh",OrientFaceSphere,Phys.RadEarth,Grid)
  elseif GridType == "CubedSphere"
    Grid = CubedGrid(nPanel,OrientFaceSphere,Phys.RadEarth,Grid)
  elseif GridType == "TriangularSphere"
    IcosahedronGrid = CGDycore.CreateIcosahedronGrid()
    RefineLevel =  0
    for iRef = 1 : RefineLevel
      CGDycore.RefineEdgeTriangularGrid!(IcosahedronGrid)
      CGDycore.RefineFaceTriangularGrid!(IcosahedronGrid)
    end
    CGDycore.NumberingTriangularGrid!(IcosahedronGrid)
    Grid = CGDycore.TriangularGridToGrid(IcosahedronGrid,Rad,Grid)
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

  if Model.Stretch
    if Model.StretchType == "ICON"  
      sigma = 1.0
      lambda = 3.16
      AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
    elseif Model.StretchType == "Exp"  
      AddExpStretchVerticalGrid!(SubGrid,nz,H,30.0,1500.0)
    else
      AddVerticalGrid!(SubGrid,nz,H)
    end
  else  
    AddVerticalGrid!(SubGrid,nz,H)
  end

  Exchange = ExchangeStruct(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Model.HorLimit)
  Output = OutputStruct(Topography)
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,SubGrid,Model,TimeStepper,ParallelCom,Output,Exchange,DoF,nz,
    Model.NumV,Model.NumTr,())
  CG = CGStruct(backend,FT,OrdPoly,OrdPolyZ,Global.Grid)
  (CG,Metric) = DiscretizationCG(backend,FT,JacobiSphere3,CG,Global)

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = vtkStruct{FT}(backend,1,TransSphereX!,CG,Metric,Global)
  unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp

  if Topography.TopoS == "EarthOrography"
    zS = Orography(CG,Global)
    Output.RadPrint = H
    Output.Flat=false
    nzTemp = Global.Grid.nz
    Global.Grid.nz = 1
    vtkCacheOrography = vtkStruct(OrdPoly,TransSphereX,CG,Global)
    unstructured_vtkOrography(zS,vtkCacheOrography, Global.Grid.NumFaces, CG,  Proc, ProcNumber)
    Global.Grid.nz = nzTemp
  end

  if Topography.TopoS == "EarthOrography"
    (CG,Metric) = DiscretizationCG(backend,FT,JacobiSphere3,CG,Global,zS)
  else
    (CG,Metric) = DiscretizationCG(backend,FT,JacobiSphere3,CG,Global)
  end
  return CG,Metric,Global
end  


function InitCart(backend,FT,OrdPoly,OrdPolyZ,nx,ny,Lx,Ly,x0,y0,nz,H,Boundary,GridType,Topography,Decomp,Model,Phys)    

  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1 
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = ParallelComStruct()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)
  Grid = GridStruct(nz,Topography)

  Grid = CartGrid(nx,ny,Lx,Ly,x0,y0,OrientFaceCart,Boundary,Grid)

  CellToProc = Decompose(Grid,ProcNumber)
  SubGrid = ConstructSubGrid(Grid,CellToProc,Proc)

  @show "Stretch"
  if Model.Stretch
    if Model.StretchType == "ICON"  
      sigma = 1.0
      lambda = 3.16
      AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
    elseif Model.StretchType == "Exp"  
      AddExpStretchVerticalGrid!(SubGrid,nz,H,30.0,1500.0)
    else
      AddVerticalGrid!(SubGrid,nz,H)
    end
  else  
    AddVerticalGrid!(SubGrid,nz,H)
  end

  Exchange = ExchangeStruct(SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Model.HorLimit)
  Output = OutputStruct(Topography)
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,SubGrid,Model,TimeStepper,ParallelCom,Output,Exchange,DoF,nz,
    Model.NumV,Model.NumTr,())
  CG = CGStruct(backend,FT,OrdPoly,OrdPolyZ,Global.Grid)
  (CG,Metric) = DiscretizationCG(backend,FT,JacobiDG3,CG,Global)

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = vtkStruct{FT}(backend,1,TransCartX!,CG,Metric,Global)
  unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp

  if Topography.TopoS == "EarthOrography"
    zS = Orography(CG,Global)
    Output.RadPrint = H
    Output.Flat=false
    nzTemp = Global.Grid.nz
    Global.Grid.nz = 1
    vtkCacheOrography = vtkStruct{FT}(backend,OrdPoly,TransCartX!,CG,Metric,Global)
    unstructured_vtkOrography(zS,vtkCacheOrography, Global.Grid.NumFaces, CG,  Proc, ProcNumber)
    Global.Grid.nz = nzTemp
    (CG,Metric) = DiscretizationCG(backend,FT,JacobiDG3,CG,Global,zS)
  else
    (CG,Metric) = DiscretizationCG(backend,FT,JacobiDG3,CG,Global)
  end



  return CG, Metric, Global
end  


