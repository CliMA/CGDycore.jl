# Load ETOPO1 ice-sheet surface data
# Ocean values are considered 0
data = NCDataset("/groups/esm/asridhar/ETOPO1_Ice_g_gdal.grd")
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
function earth_orography(coords)
  λ, ϕ = coords.long, coords.lat # Unpack longitude
  FT = eltype(λ)
  zₛ = spline_2d(λ,ϕ)
  return FT(zₛ)
end
