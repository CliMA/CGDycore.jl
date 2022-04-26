#function testStaticArrays

using CGDycore
using StaticArrays

OrdPoly = 4
OrdPolyZ=1
nz = 20
nPanel = 16
NF = 6 * nPanel * nPanel
NumV = 5
const NPOLY = 5
const NZ = 20

Field = Vector{CGDycore.ColumnElementStruct}(undef, 5)
#Field = Vector{CGDycore.ColumnElementStruct}(CGDycore.ColumnElement(NPOLY,NZ,Float64), 5)
for i=1:5
  Field[i] = CGDycore.ColumnElement(NPOLY,NZ,Float64)
end  
Column=CGDycore.ColumnElement(NPOLY,NZ,Float64)
Column.Rho .= 1.0
Column.Rho[2,7] = 5.0

@SArray b=zeros(Float64,NPOLY,NPOLY,NZ)
@time @views CGDycore.test1!(Field[1],Column,NPOLY,NZ)


# Physical parameters

