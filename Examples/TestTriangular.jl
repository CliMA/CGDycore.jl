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
CGDycore.NumberingTriangularGrid!(IcosahedronGrid)
nz = 1
Rad = 1.0
Topography = (TopoS=TopoS,H=1)

GridTri = CGDycore.GridStruct(nz,Topography)
GridTri = CGDycore.TriangularGridToGrid(IcosahedronGrid,Rad,GridTri)
vtkSkeletonMeshTri = CGDycore.vtkStruct{Float64}(backend,GridTri)
CGDycore.vtkSkeleton(vtkSkeletonMeshTri,"IcosahedronTri", 1, 1)

GridDel = CGDycore.GridStruct(nz,Topography)
GridDel = CGDycore.DelaunayGridToPolyGrid(IcosahedronGrid,Rad,GridDel)
vtkSkeletonMeshDel = CGDycore.vtkStruct{Float64}(backend,GridDel)
CGDycore.vtkSkeleton(vtkSkeletonMeshDel,"IcosahedronDel", 1, 1)


# RefineCellTriangularGrid(IcosahedronGrid)

# CALL NumberingTriangularGrid(Square)
# CALL TriangularGridToPolyGrid(nz,Square,PolyGrid)

a = 4

