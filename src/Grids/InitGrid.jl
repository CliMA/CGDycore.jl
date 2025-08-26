function InitGridSphere(backend,FT,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,GridType,Decomp,RadEarth,Model,
  ParallelCom;order=true,ChangeOrient=3,Discretization="CG")

  ProcNumber = ParallelCom.ProcNumber
  Proc = ParallelCom.Proc

  if GridType == "HealPix"
  # Grid=CGDycore.InputGridH("Grid/mesh_H12_no_pp.nc",
  # CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)
  # Grid=Grids.InputGridH(backend,FT,"Grid/mesh_H24_no_pp.nc", Grids.OrientFaceSphere,RadEarth,nz)
    Grid = Grids.HealpixGrid(backend,FT,ns,RadEarth,nz)
  elseif GridType == "SQuadGen"
    Grid = Grids.InputGrid(backend,FT,"Grid/baroclinic_wave_2deg_x4.g",Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "Msh"
#   Grid = Grids.InputGridMsh(backend,FT,"Grid/natural_earth.msh",Grids.OrientFaceSphere,RadEarth,nz)
    Grid = Grids.InputGridMsh(backend,FT,"Grid/Quad.msh",Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "CubedSphere"
    Grid = Grids.CubedGrid(backend,FT,nPanel,Grids.OrientFaceSphere,RadEarth,nz,order=order)
  elseif GridType == "TriangularSphere"
    Grid = TriangularGrid(backend,FT,RefineLevel,RadEarth,nz;ChangeOrient=ChangeOrient)
  elseif GridType == "DelaunaySphere"
    Grid = DelaunayGrid(backend,FT,RefineLevel,RadEarth,nz)
  elseif GridType == "MPASO"
    Grid=Grids.InputGridMPASO(backend,FT,"Grid/QU240.nc", Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "MPAS"
    Grid=Grids.InputGridMPASO(backend,FT,"Grid/x4.163842.grid.nc", Grids.OrientFaceSphere,RadEarth,nz)
#   Grid=Grids.InputGridMPASO(backend,FT,"Grid/x1.40962.grid.nc", Grids.OrientFaceSphere,RadEarth,nz)
  elseif GridType == "MPAS2Triangle"
    Grid=Grids.InputGridMPAS2Triangular(backend,FT,"Grid/x1.2562.grid.nc", 
      Grids.OrientFaceSphere,RadEarth,nz;ChangeOrient)
  elseif GridType == "SphericalGrid"
    Grid = SphericalGrid(backend,FT,nLon,nLat,LatB,Grids.OrientFaceSphere,RadEarth,nz;ChangeOrient)
  elseif GridType == "ICON"
    Grid=Grids.InputGridICON(backend,FT,"Grid/icon_grid_0009_R02B03_R.nc", Grids.OrientFaceSphere,
      RadEarth,nz;ChangeOrient)
  else
    @show "False GridType"
  end
  if Proc == 1
    @show Grid.NumFaces   
    @show Grid.NumEdges   
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
  SubGrid = Grids.ConstructSubGridGhost(Grid,CellToProc,Proc,order=order)

#=
  if Discretization == "DG"
    Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly+1,0,CellToProc,Proc,ProcNumber,
      Model.HorLimit;Discretization="DG")
  elseif Discretization == "CG"
    Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly-1,1,CellToProc,Proc,ProcNumber,
    Model.HorLimit)
  end
=#  
  return SubGrid, CellToProc
end  

function InitGridCart(backend,FT,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom;
  order=true,GridType="Quad",Discretization="CG",ChangeOrient=3)

  ProcNumber = ParallelCom.ProcNumber
  Proc = ParallelCom.Proc

  if GridType == "Quad"
    Grid = Grids.CartGrid(backend,FT,nx,ny,Lx,Ly,x0,y0,Grids.OrientFaceCart,Boundary,nz;order)
  elseif GridType == "Tri"  
    Grid = Grids.CartGridTri(backend,FT,nx,ny,Lx,Ly,x0,y0,Grids.OrientFaceCart,Boundary,nz;order,ChangeOrient=ChangeOrient)
  else  
    stop  
  end  
  CellToProc = Grids.Decompose(Grid,nx,ny,ProcNumber)
  SubGrid = Grids.ConstructSubGridGhost(Grid,CellToProc,Proc;order)
#=  
  if Discretization == "DG"
    Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly+1,0,CellToProc,Proc,ProcNumber,Model.HorLimit;
      Discretization="DG")
  else  
    Exchange = Parallels.ExchangeStruct{FT}(backend,SubGrid,OrdPoly-1,1,CellToProc,Proc,ProcNumber,Model.HorLimit)
  end  
=#
  return SubGrid, CellToProc
end  
