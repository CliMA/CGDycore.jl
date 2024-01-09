function InitSphere(backend,FT,OrdPoly,OrdPolyZ,nz,nPanel,H,GridType,Topography,Decomp,Model,Phys,RadEarth)    

  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1 
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = ParallelComStruct()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)
# Grid = Grids.GridStruct{FT}(backend,nz,Topography)

  if GridType == "HealPix"
  # Grid=CGDycore.InputGridH("Grid/mesh_H12_no_pp.nc",
  # CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
    Grid=Grids.InputGridH(backend,FT,"Grid/mesh_H24_no_pp.nc", Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "SQuadGen"
    Grid = Grids.InputGrid("Grid/baroclinic_wave_2deg_x4.g",Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "Msh"
    Grid = Grids.InputGridMsh(backend,FT,"Grid/Quad.msh",Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "CubedSphere"
    Grid = Grids.CubedGrid(backend,FT,nPanel,Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "TriangularSphere"
    IcosahedronGrid = Grids.CreateIcosahedronGrid()
    RefineLevel =  0
    for iRef = 1 : RefineLevel
      Grids.RefineEdgeTriangularGrid!(IcosahedronGrid)
      Grids.RefineFaceTriangularGrid!(IcosahedronGrid)
    end
    Grids.NumberingTriangularGrid!(IcosahedronGrid)
    Grid = Grids.TriangularGridToGrid(IcosahedronGrid,RadEarth,Grid)
  else
    @show "False GridType"
  end

  if Decomp == "Hilbert"
    Parallels.HilbertFaceSphere!(Grid)
    CellToProc = Grids.Decompose(Grid,ProcNumber)
  elseif Decomp == "EqualArea"
    CellToProc = Grids.DecomposeEqualArea(Grid,ProcNumber)
  else
    CellToProc = ones(Int,Grid.NumFaces)
    println(" False Decomp method ")
  end

  SubGrid = Grids.ConstructSubGrid(Grid,CellToProc,Proc)

  if Model.Stretch
    if Model.StretchType == "ICON"  
      sigma = 1.0
      lambda = 3.16
      Grids.AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
    elseif Model.StretchType == "Exp"  
      Grids.AddExpStretchVerticalGrid!(SubGrid,nz,H,30.0,1500.0)
    else
      Grids.AddVerticalGrid!(SubGrid,nz,H)
    end
  else  
    Grids.AddVerticalGrid!(SubGrid,nz,H)
  end

  Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Model.HorLimit)
  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,SubGrid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr,())
  CG = CGQuad{FT}(backend,OrdPoly,OrdPolyZ,Global.Grid)

  if Topography.TopoS == "EarthOrography"
    zS = Grids.Orography(backend,FT,CG,Exchange,Global)
  end

  if Topography.TopoS == "EarthOrography"
    (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiSphere3,CG,Exchange,Global,zS)
  else
    (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiSphere3,CG,Exchange,Global)
  end
  Global.Output.dTol = 2*pi / nPanel
  # Output Orography
  if Topography.TopoS == "EarthOrography"
    Output.Flat=true
    nzTemp = Global.Grid.nz
    Global.Grid.nz = 1
    vtkCacheOrography = Outputs.vtkInit2D(CG.OrdPoly,Grids.TransSphereX!,CG,Metric,Global)
    Outputs.unstructured_vtkOrography(zS,vtkCacheOrography,Global.Grid.NumFaces,CG,Proc,ProcNumber)
    Global.Grid.nz = nzTemp
  end
  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,1,Grids.TransSphereX!,CG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp
  return CG,Metric,Exchange,Global
end  


function InitCart(backend,FT,OrdPoly,OrdPolyZ,nx,ny,Lx,Ly,x0,y0,nz,H,Boundary,GridType,Decomp,Model,Phys,TopoProfile)

  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1 
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = ParallelComStruct()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Grid = Grids.CartGrid(backend,FT,nx,ny,Lx,Ly,x0,y0,Grids.OrientFaceCart,Boundary,nz)

  CellToProc = Grids.Decompose(Grid,ProcNumber)
  SubGrid = Grids.ConstructSubGrid(Grid,CellToProc,Proc)

  if Model.Stretch
    if Model.StretchType == "ICON"  
      sigma = 1.0
      lambda = 3.16
      Grids.AddStretchICONVerticalGrid!(SubGrid,nz,H,sigma,lambda)
    elseif Model.StretchType == "Exp"  
      Grids.AddExpStretchVerticalGrid!(SubGrid,nz,H,30.0,1500.0)
    else
      Grids.AddVerticalGrid!(SubGrid,nz,H)
    end
  else  
    Grids.AddVerticalGrid!(SubGrid,nz,H)
  end

  Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Model.HorLimit)
  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,SubGrid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr,())
  CG = CGQuad{FT}(backend,OrdPoly,OrdPolyZ,Global.Grid)
  (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiDG3,CG,Exchange,Global)

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,1,Grids.TransCartX!,CG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart, Global.Grid.NumFaces, Proc, ProcNumber)
  Global.Grid.nz = nzTemp

  zS = Grids.Orography(backend,FT,CG,Global,TopoProfile)
  @show minimum(zS),maximum(zS)
  (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiDG3,CG,Exchange,Global,zS)

  return CG, Metric, Exchange, Global
end  


