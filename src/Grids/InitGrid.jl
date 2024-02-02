function InitGridSphere(backend,FT,OrdPoly,nz,nPanel,RefineLevel,GridType,Decomp,RadEarth,Model,ParallelCom)    

  ProcNumber = ParallelCom.ProcNumber
  Proc = ParallelCom.Proc

  if GridType == "HealPix"
  # Grid=CGDycore.InputGridH("Grid/mesh_H12_no_pp.nc",
  # CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
    Grid=Grids.InputGridH(backend,FT,"Grid/mesh_H24_no_pp.nc", Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "SQuadGen"
    Grid = Grids.InputGrid(backend,FT,"Grid/baroclinic_wave_2deg_x4.g",Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "Msh"
    Grid = Grids.InputGridMsh(backend,FT,"Grid/natural_earthNeu.msh",Grids.OrientFaceSphere,RadEarth,nz)
#   Grid = Grids.InputGridMsh(backend,FT,"Grid/Quad.msh",Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "CubedSphere"
    Grid = Grids.CubedGrid(backend,FT,nPanel,Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "TriangularSphere"
    Grid = TriangularGrid(backend,FT,RefineLevel,RadEarth,nz)
  elseif GridType == "DelaunaySphere"
    Grid = DelaunayGrid(backend,FT,RefineLevel,RadEarth,nz)
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

  Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly,CellToProc,Proc,ProcNumber,Model.HorLimit)
  return SubGrid, Exchange
end  
