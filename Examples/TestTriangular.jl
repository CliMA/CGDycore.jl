#import src/Grid/Triangular.jl
#include("../src/Grid/Triangular.jl")
include("../src/CGDycore.jl")
using PairedLinkedLists
using KernelAbstractions

backend = CPU()

parsed_args = CGDycore.parse_commandline()
TopoS = parsed_args["TopoS"]

IcosahedronGrid = CGDycore.CreateIcosahedronGrid()
RefineLevel =  3
for iRef = 1 : RefineLevel
  CGDycore.RefineEdgeTriangularGrid!(IcosahedronGrid)  
  CGDycore.RefineFaceTriangularGrid!(IcosahedronGrid)
end
nz = 1
Topography = (TopoS=TopoS,H=1)
Grid = CGDycore.GridStruct(nz,Topography)
Rad = 1.0
Grid = CGDycore.TriangularGridToGrid(IcosahedronGrid,Rad,Grid)
vtkSkeletonMesh = CGDycore.vtkStruct{Float64}(backend,Grid)
CGDycore.vtkSkeleton(vtkSkeletonMesh,"Icosahedron", 1, 1)


# RefineCellTriangularGrid(IcosahedronGrid)

# CALL NumberingTriangularGrid(Square)
# CALL TriangularGridToPolyGrid(nz,Square,PolyGrid)

a = 4

