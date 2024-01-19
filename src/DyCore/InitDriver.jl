function InitSphere(backend,FT,OrdPoly,OrdPolyZ,H,Topography,Model,Phys,TopoProfile,Exchange,Grid,ParallelCom)    
  nz = Grid.nz 
  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr,())
  CG = CGQuad{FT}(backend,OrdPoly,OrdPolyZ,Global.Grid)

  if Model.Stretch
    if Model.StretchType == "ICON"
      sigma = 1.0
      lambda = 3.16
      Grids.AddStretchICONVerticalGrid!(Global.Grid,nz,H,sigma,lambda)
    elseif Model.StretchType == "Exp"
      Grids.AddExpStretchVerticalGrid!(Global.Grid,nz,H,30.0,1500.0)
    else
      Grids.AddVerticalGrid!(Global.Grid,nz,H)
    end
  else
    Grids.AddVerticalGrid!(Global.Grid,nz,H)
  end

  if Topography.TopoS == "EarthOrography"
    zS = Grids.Orography(backend,FT,CG,Exchange,Global)
  else
    zS = Grids.Orography(backend,FT,CG,Exchange,Global,TopoProfile)
  end


  (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiSphere3,CG,Exchange,Global,zS)

  # Output Orography
  Global.Output.dTol = 2*pi / 30
# Global.Output.dTol = 1.e-8
  Output.Flat=true
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCacheOrography = Outputs.vtkInit2D(CG.OrdPoly,Grids.TransSphereX!,CG,Metric,Global)
  Outputs.unstructured_vtkOrography(zS,vtkCacheOrography,Global.Grid.NumFaces,CG,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,1,Grids.TransSphereX!,CG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  return CG, Metric, Global
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

  zS = Grids.Orography(backend,FT,CG,Exchange,Global,TopoProfile)
  (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiDG3,CG,Exchange,Global,zS)

  return CG, Metric, Exchange, Global
end  


