using NCDatasets
using DelimitedFiles
using Dierckx
using CGDycore

OrdPoly = 4
nz = 40

OrdPolyZ=1
nPanel = 8
NF = 6 * nPanel * nPanel
NumV = 5
NumTr = 2

# Physical parameters
Phys=CGDycore.PhysParameters()

# Grid
H = 30000.0
#H = 45000.0
Topography=(TopoS="",H=H,Rad=Phys.RadEarth)
Grid=CGDycore.Grid(nz,Topography)
Grid=CGDycore.CubedGrid(nPanel,CGDycore.OrientFaceSphere,Phys.RadEarth,Grid)

(w,xw)=CGDycore.GaussLobattoQuad(OrdPoly)
lam = zeros(OrdPoly+1,OrdPoly+1,NF)
phi = zeros(OrdPoly+1,OrdPoly+1,NF)
for iF = 1 : NF
  F = Grid.Faces[iF]  
  for j = 1 : OrdPoly+1
    for i = 1 : OrdPoly+1
      X1=0.25*((1-xw[i])*(1-xw[j])*F.P[1].x+
       (1+xw[i])*(1-xw[j])*F.P[2].x+
       (1+xw[i])*(1+xw[j])*F.P[3].x+
       (1-xw[i])*(1+xw[j])*F.P[4].x);
      X2=0.25*((1-xw[i])*(1-xw[j])*F.P[1].y+
       (1+xw[i])*(1-xw[j])*F.P[2].y+
       (1+xw[i])*(1+xw[j])*F.P[3].y+
       (1-xw[i])*(1+xw[j])*F.P[4].y);
      X3=0.25*((1-xw[i])*(1-xw[j])*F.P[1].z+
       (1+xw[i])*(1-xw[j])*F.P[2].z+
       (1+xw[i])*(1+xw[j])*F.P[3].z+
       (1-xw[i])*(1+xw[j])*F.P[4].z);
      (lam[i,j,iF], phi[i,j,iF])  = CGDycore.cart2sphere(X1,X2,X3)
    end
  end
end  


# Load ETOPO1 ice-sheet surface data
# Ocean values are considered 0
data = NCDataset("/Users/knoth/Documents/GitHub/CliMA/CGDycore.jl/ETOPO1_Ice_g_gdal.grd")
# Unpack information
x_range = data["x_range"][:] 
y_range = data["y_range"][:]
z_range = data["z_range"][:]
spacing = data["spacing"][:]
dimension = data["dimension"][:]
elevation = data["z"][:]
lon = collect(x_range[1]:spacing[1]:x_range[2])
lat = collect(y_range[1]:spacing[2]:y_range[2])
nlon = dimension[1]
nlat = dimension[2]
zlevels = reshape(elevation, (nlon, nlat))

# Construct Elevation Matrix
function coarsen(X::Vector; factor = 100)
  return X[1:factor:end]
end
function coarsen(X::AbstractArray; factor = 100)
  return X[1:factor:end, 1:factor:end]
end
map_source = zlevels
map_source[map_source .< 0.0] .= 0.0
#const spline_2d=Spline2D(coarsen(lon), coarsen(lat), reverse(coarsen(map_source), dims=2);kx=3, ky=3, s=0.0)
const spline_2d=Spline2D(lon, lat, reverse(map_source, dims=2);kx=3, ky=3, s=0.0)

"""
  earth_orography(λ,Φ)
λ = Longitude (degrees)
ϕ = Latitude (degrees)
Simple spline interpolant representation of Earth's surface
orography. Data sourced from ETOPO1 datasets (grid-referenced),
sea-ice surface information.
"""
function earth_orography(lambda,phi)
  zs = evaluate(spline_2d,lambda,phi)
end

zs = zeros(OrdPoly+1,OrdPoly+1,NF)
for iF = 1 : NF
  for j = 1 : OrdPoly+1
    for i = 1 : OrdPoly+1
      lamD = lam[i,j,iF] * 180.0 / pi  
      phiD = phi[i,j,iF] * 180.0 / pi
      zs[i,j,iF] = earth_orography(lamD,phiD);
    end
  end
end  
