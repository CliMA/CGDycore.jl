#function testDecompose

using CGDycore

OrdPoly = 4
nz = 3

OrdPolyZ=1
nPanel = 2
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 2


# Physical parameters
Phys=CGDycore.PhysParameters()

Param=(T0E=310.0,
      )
Model = CGDycore.Model(Param)

# Grid
H = 30000.0
#H = 45000.0
Topography=(TopoS="EarthOrography",H=H,Rad=Phys.RadEarth)
Output=CGDycore.Output(Topography)

ProcNumber = 2

Proc = 1
Grid1=CGDycore.Grid(nz,Topography)
Grid1=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid1)
CellToProc = CGDycore.Decompose(Grid1,ProcNumber)
SubGrid1 = CGDycore.ConstructSubGrid(Grid1,CellToProc,Proc)
CGDycore.AddVerticalGrid!(SubGrid1,nz,H)
Global1 = CGDycore.Global(SubGrid1,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())
Global1.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid1.NumFaces,nz)
(CG1,Global1)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global1)
@show CG1.NumG
@show size(CG1.Boundary)

Proc = 2
Grid2=CGDycore.Grid(nz,Topography)
Grid2=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid2)
CellToProc = CGDycore.Decompose(Grid2,ProcNumber)
SubGrid2 = CGDycore.ConstructSubGrid(Grid2,CellToProc,Proc)
CGDycore.AddVerticalGrid!(SubGrid2,nz,H)
Global2 = CGDycore.Global(SubGrid2,Model,Phys,Output,OrdPoly+1,nz,NumV,NumTr,())
Global2.Metric=CGDycore.Metric(OrdPoly+1,OrdPolyZ+1,SubGrid2.NumFaces,nz)
(CG2,Global2)=CGDycore.Discretization(OrdPoly,OrdPolyZ,CGDycore.JacobiSphere3,Global2)
@show CG2.NumG
@show size(CG2.Boundary)

a = 5
