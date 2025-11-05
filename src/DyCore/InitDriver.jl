function InitSphereDG2(backend,FT,OrdPoly,OrdPolyZ,H,Topography,Model,Phys,TopoProfile,Exchange,Grid,ParallelCom)    
  nz = Grid.nz 
  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr)
  DG = FiniteElements.DGQuad{FT}(backend,OrdPoly,OrdPolyZ,Global.Grid,ParallelCom.Proc)

  (DG,Metric) = DiscretizationDG2(backend,FT,Grids.JacobiSphereDG2GPU!,DG,Exchange,Global)

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,1,Grids.TransSphereX!,DG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  return DG, Metric, Global
end  

function InitSphereDG(backend,FT,OrdPoly,OrdPolyZ,DGMethod,OrdPrint,OrdPrintZ,H,Topography,Model,
  Phys,TopoProfile,CellToProc,Grid,ParallelCom)
  nz = Grid.nz
  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr)
  if Grid.Type == Grids.Quad()
    DG = FiniteElements.DGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Global.Grid,ParallelCom.Proc)
  else
    DG = FiniteElements.DGTri{FT}(backend,DGMethod,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)
  end    
  Exchange = Parallels.ExchangeStruct{FT}(backend,Grid,DG,CellToProc,Proc,ProcNumber,
      Model.HorLimit;Discretization="DG")

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
    @time zS, GradDx_zs, GradDy_zs = Grids.Orography4(backend,FT,CG,Exchange,Global)
  else
    zS = Grids.Orography(backend,FT,DG,Exchange,Global,TopoProfile,Grid.Type)
  end

  Metric = DiscretizationDG(backend,FT,Grids.JacobiSphereDG3GPU!,DG,Exchange,Global,zS,Grid.Type)

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,0,0,Grids.TransSphereX!,DG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  return DG, Metric, Exchange, Global
end

function InitSphere(backend,FT,OrdPoly,OrdPolyZ,OrdPrint,H,Topography,Model,Phys,TopoProfile,
  CellToProc,Grid,ParallelCom)
  nz = Grid.nz 
  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr)
  CG = FiniteElements.CGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,Global.Grid)
  Exchange = Parallels.ExchangeStruct{FT}(backend,Grid,CG,CellToProc,Proc,ProcNumber,
      Model.HorLimit;Discretization="CG")

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
    #zS, GradDx_zs, GradDy_zs = Grids.Orography4(backend,FT,CG,Exchange,Global)
#   zS = Grids.Orography2(backend,FT,CG,Exchange,Global)
    zS = Grids.Orography4(backend,FT,CG,Exchange,Global)
  else
    zS = Grids.Orography(backend,FT,CG,Exchange,Global,TopoProfile,Grid.Type)
  end


  (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiSphere3,CG,Exchange,Global,zS)

  # Output Orography
  Global.Output.dTol = 2*pi / 30
# Global.Output.dTol = 1.e-8
# Output.Flat=true
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCacheOrography = Outputs.vtkInit2D(CG.OrdPoly,Grids.TransSphere2DX!,CG,Metric,Global)
  Outputs.unstructured_vtkOrography(zS,vtkCacheOrography,Global.Grid.NumFaces,CG,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,0,0,Grids.TransSphereX!,CG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  return CG, Metric, Exchange, Global
end  

function InitCart(backend,FT,OrdPoly,OrdPolyZ,OrdPrint,H,Topography,Model,Phys,TopoProfile,
  CellToProc,Grid,ParallelCom)

  nz = Grid.nz
  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr)
  CG = FiniteElements.CGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,Global.Grid)
  Exchange = Parallels.ExchangeStruct{FT}(backend,Grid,CG,CellToProc,Proc,ProcNumber,
      Model.HorLimit;Discretization="CG")

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
    zS = Grids.Orography(backend,FT,CG,Exchange,Global,TopoProfile,Grid.Type)
  end

  (CG,Metric) = DiscretizationCG(backend,FT,Grids.JacobiDG3,CG,Exchange,Global,zS)

  # Output Orography
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCacheOrography = Outputs.vtkInit2D(CG.OrdPoly,Grids.TransCartX!,CG,Metric,Global)
  Outputs.unstructured_vtkOrography(zS,vtkCacheOrography,Global.Grid.NumFaces,CG,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,0,0,Grids.TransCartX!,CG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp
  nzTemp = Global.Grid.nz

  return CG, Metric, Exchange, Global
end  

function InitCartDG(backend,FT,OrdPoly,OrdPolyZ,DGMethod,OrdPrint,OrdPrintZ,H,Topography,Model,
  Phys,TopoProfile,CellToProc,Grid,ParallelCom)
  nz = Grid.nz
  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber

  TimeStepper = TimeStepperStruct{FT}(backend)

  Output = OutputStruct()
  DoF = (OrdPoly + 1) * (OrdPoly + 1)
  Global = GlobalStruct{FT}(backend,Grid,Model,TimeStepper,ParallelCom,Output,DoF,nz,
    Model.NumV,Model.NumTr)
  if Grid.Type == Grids.Quad()
    DG = FiniteElements.DGQuad{FT}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Global.Grid,ParallelCom.Proc)
  else
    DG = FiniteElements.DGTri{FT}(backend,DGMethod,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)
  end    
  Exchange = Parallels.ExchangeStruct{FT}(backend,Grid,DG,CellToProc,Proc,ProcNumber,
      Model.HorLimit;Discretization="DG")

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
    zS = Grids.Orography(backend,FT,DG,Exchange,Global,TopoProfile,Grid.Type)
  end

  Metric = DiscretizationDG(backend,FT,Grids.JacobiCartDG3GPU!,DG,Exchange,Global,zS,Grid.Type)

  # Output Orography
# nzTemp = Global.Grid.nz
# Global.Grid.nz = 1
# vtkCacheOrography = Outputs.vtkInit2D(DG.OrdPoly,Outputs.TransCartX!,DG,Metric,Global)
# Outputs.unstructured_vtkOrography(zS,vtkCacheOrography,Global.Grid.NumFaces,DG,Proc,ProcNumber)
# Global.Grid.nz = nzTemp

  # Output partition
  nzTemp = Global.Grid.nz
  Global.Grid.nz = 1
  vtkCachePart = Outputs.vtkStruct{FT}(backend,0,0,Grids.TransCartX!,DG,Metric,Global)
  Outputs.unstructured_vtkPartition(vtkCachePart,Global.Grid.NumFaces,Proc,ProcNumber)
  Global.Grid.nz = nzTemp

  return DG, Metric, Exchange, Global
end  


