using NCDatasets
using StructArrays

mutable struct TilesRawGrid_T
  altitude::Array{Union{Missing, Int32}, 2}
  lon::Array{Float32, 1}
  lat::Array{Float32, 1}
  lon_max::Float64
  lon_min::Float64
  lat_max::Float64
  lat_min::Float64
  start_lon::Float64
  end_lon::Float64
  start_lat::Float64
  end_lat::Float64
  dlon::Float64
  dlat::Float64
  nlon::Int
  nlat::Int
end

function TilesRawGrid_T()
  altitude = zeros(Int32,0,0)
  lon = zeros(Float32,0)
  lat = zeros(Float32,0)
  lon_max = 0.0
  lon_min = 0.0
  lat_max = 0.0
  lat_min = 0.0
  start_lon = 0.0
  end_lon = 0.0
  start_lat = 0.0
  end_lat = 0.0
  dlon = 0.0
  dlat = 0.0
  nlon = 0
  nlat = 0
  return TilesRawGrid_T(
    altitude,
    lon,
    lat,
    lon_max,
    lon_min,
    lat_max,
    lat_min,
    start_lon,
    end_lon,
    start_lat,
    end_lat,
    dlon,
    dlat,
    nlon,
    nlat)
end  

function TopoData()
  nlon = 200
  nlat = 100
  dlon = 360 / nlon
  dlat = 180 / nlat
  lon = collect(0 : dlon : 360)
  lat = collect(90 : dlat : -90)
  nlon = nlon + 1
  nlat = nlat + 1
  nlon1 = 1
  nlon2 = nlon
  nlat1 = 51-0
  nlat2 = 51+0
  zlevels = zeros(Float64,nlon,nlat)
  for ilat = nlat1 : nlat2
    zlevels[nlon1:nlon2,ilat] .= 200.0 
  end  
  return (lon, lat, zlevels)
end
function TopoDataETOPO(MinLonL,MaxLonL,MinLonR,MaxLonR,MinLat,MaxLat)
  # Load ETOPO1 ice-sheet surface data
  # Ocean values are considered 0
  ds = NCDataset("/Users/knoth/Documents/GitHub/CliMA/CGDycore.jl/ETOPO1_Ice_g_gdal.grd")
  # Unpack information
  x_range = ds["x_range"][:]
  y_range = ds["y_range"][:]
  z_range = ds["z_range"][:]
  spacing = ds["spacing"][:]
  dimension = ds["dimension"][:]
  elevation = ds["z"][:]
  lon = collect(x_range[1]:spacing[1]:x_range[2])
  lat = collect(y_range[1]:spacing[2]:y_range[2])
  nlon = dimension[1]
  nlat = dimension[2]
  dLon = 360.0 / nlon
  dLat = 180.0 / nlat
  ilonLS = max(floor(Int,(MinLonL+180.)/dLon),1)
  ilonLE = min(ceil(Int,(MaxLonL+180.)/dLon),nlon)
  ilonRS = max(floor(Int,(MinLonR+180.)/dLon),1)
  ilonRE = min(ceil(Int,(MaxLonR+180.)/dLon),nlon)
  ilatS = max(floor(Int,(MinLat+90.)/dLat),1)
  ilatE = min(ceil(Int,(MaxLat+90.)/dLat),nlat)
  zlevels = max.(reshape(elevation, (nlon, nlat)), 0.0)
  return (lon[ilonLS:ilonLE], lon[ilonRS:ilonRE], lat[ilatS:ilatE], 
    zlevels[ilonLS:ilonLE,ilatS:ilatE], zlevels[ilonRS:ilonRE,ilatS:ilatE])
end  

function TopoDataGLOBE()
  deg2rad = pi / 180.0
  RadEarth = 6.37122e+6
  ntiles_column = 4
  ntiles_row = 4
  ntiles = ntiles_column * ntiles_row
  list=["GLOBE_A10.nc",
        "GLOBE_B10.nc",
        "GLOBE_C10.nc",
        "GLOBE_D10.nc",
        "GLOBE_E10.nc",
        "GLOBE_F10.nc",
        "GLOBE_G10.nc",
        "GLOBE_H10.nc",
        "GLOBE_I10.nc",
        "GLOBE_J10.nc",
        "GLOBE_K10.nc",
        "GLOBE_L10.nc",
        "GLOBE_M10.nc",
        "GLOBE_N10.nc",
        "GLOBE_O10.nc",
        "GLOBE_P10.nc"]

  TilesRawGrid = map(1:ntiles) do i
    TilesRawGrid_T()
  end
  i = 1
  nlon_tot = 0
  nlat_tot = 0
  for file in list
    ds = NCDataset("/Users/knoth/Documents/GitHub/CliMA/CGDycore.jl/Topo/"*file)  
    nlon = ds.dim["lon"]
    nlat = ds.dim["lat"]
    nlon_tot += nlon
    nlat_tot += nlat
    close(ds)
#   TilesRawGrid[i].altitude = data["altitude"][:,:]
#   TilesRawGrid[i].altitude.attrib["_FillValue"] = 0
#   TilesRawGrid[i].lon = data["lon"][:]
#   TilesRawGrid[i].lat = data["lat"][:]
#   TilesRawGrid[i].lon_min = minimum(TilesRawGrid[i].lon)
#   TilesRawGrid[i].lon_max = maximum(TilesRawGrid[i].lon)
#   TilesRawGrid[i].lat_min = minimum(TilesRawGrid[i].lat)
#   TilesRawGrid[i].lat_max = maximum(TilesRawGrid[i].lat)
#   TilesRawGrid[i].nlon = size(TilesRawGrid[i].lon,1)
#   TilesRawGrid[i].nlat = size(TilesRawGrid[i].lat,1)
#   @show TilesRawGrid[i].nlon,TilesRawGrid[i].nlat
#   @show TilesRawGrid[i].lon_min,TilesRawGrid[i].lon_max
#   @show TilesRawGrid[i].lat_min,TilesRawGrid[i].lat_max
#   nlon_tot += TilesRawGrid[i].nlon
#   nlat_tot += TilesRawGrid[i].nlat
#   TilesRawGrid[i].dlon = (TilesRawGrid[i].lon_max - TilesRawGrid[i].lon_min) / TilesRawGrid[i].nlon
#   TilesRawGrid[i].dlat = (TilesRawGrid[i].lat_max - TilesRawGrid[i].lat_min) / TilesRawGrid[i].nlat
#   # latitude from north to south, negative increment
#   TilesRawGrid[i].start_lon  += 0.5 * TilesRawGrid[i].dlon
#   TilesRawGrid[i].end_lon  -= 0.5 * TilesRawGrid[i].dlon
#   # latitude from north to south, note the negative increment!
#   TilesRawGrid[i].start_lat  += 0.5 * TilesRawGrid[i].dlat
#   # latitude from north to south, note the negative increment!
#   TilesRawGrid[i].end_lat -= 0.5 * TilesRawGrid[i].dlat
#   @show minimum(TilesRawGrid[i].lon),maximum(TilesRawGrid[i].lon)
#   @show minimum(TilesRawGrid[i].lat),maximum(TilesRawGrid[i].lat)
#   @show size(TilesRawGrid[i].lon),size(TilesRawGrid[i].lat)
#   @show size(TilesRawGrid[i].altitude)
    println("")
#   i += 1
  end   
  @show nlon_tot
  @show nlat_tot
# Altitude = zeros(Int32,nlon_tot,nlat_tot)
    ds = NCDataset("/Users/knoth/Documents/GitHub/CliMA/CGDycore.jl/Topo/"*list[13])
    a = coalesce.(ds["altitude"],Int32(0))
    nlon = ds.dim["lon"]
    nlat = ds.dim["lat"]
    nlonA = 1
    nlonE = nlon
    nlatA = 1
    nlatE = nlat
    @show maximum(a),minimum(a)
    close(ds)
  stop
  nlonA = 0 
  nlon = TilesRawGrid[13].nlon
  nlatA = 0 
  nlat = TilesRawGrid[13].nlat
  for i = 1 : nlon
    ilon = i + nlonA  
    for j = 1 : nlat  
      jlat = j + nlatA  
      @show i,j,TilesRawGrid[13].altitude[i,j]
      Altitude[ilon,jlat] = Int(TilesRawGrid[13].altitude[i,j])
    end
  end  
end

function Orography(OrdPoly,Grid,Global)
  (MinLonL,MaxLonL,MinLonR,MaxLonR,MinLat,MaxLat) = BoundingBox(Grid)
  RadEarth = Grid.Rad
  NF = Grid.NumFaces
  OP = OrdPoly + 1
  HeightCG = zeros(Float64,OP,OP,NF)
  (lonL, lonR, lat, zLevelL, zLevelR) = TopoDataETOPO(MinLonL,MaxLonL,MinLonR,MaxLonR,MinLat,MaxLat)
# (lon, lat, zLevel) = CGDycore.TopoData()
  start_Face = 1
  (Glob,NumG) = NumberingFemCG(Grid,OrdPoly);
  Height = zeros(Float64,NumG)
  NumHeight = zeros(Float64,NumG)
  (w,xw) = GaussLobattoQuad(OrdPoly)
# LenLat = length(lat)
# LenLon = length(lon)
# dLon = 360.0 / LenLon
# dLat = 180.0 / LenLat
# ilonLS = max(floor(Int,(MinLonL+180.)/dLon),1)
# ilonLE = min(ceil(Int,(MaxLonL+180.)/dLon),LenLon)
# ilonRS = max(floor(Int,(MinLonR+180.)/dLon),1)
# ilonRE = min(ceil(Int,(MaxLonR+180.)/dLon),LenLon)
# ilatS = max(floor(Int,(MinLat+90.)/dLat),1)
# ilatE = min(ceil(Int,(MaxLat+90.)/dLat),LenLat)
  for ilat = 1 : length(lat)
    for ilon = 1 : length(lonL)
      P = Point(sphereDeg2cart(lonL[ilon],-lat[ilat],RadEarth))
      (Face_id, iPosFace_id, jPosFace_id) = walk_to_nc(P,start_Face,xw,TransSphereS,RadEarth,Grid)
      start_Face = Face_id
      Inside = InsideFace(P,Grid.Faces[start_Face],Grid)
      if Inside
        iG = Glob[iPosFace_id,jPosFace_id,Face_id]
        Height[iG] += zLevelL[ilon,ilat]
        NumHeight[iG] += 1
      end  
    end
    for ilon = 1 : length(lonR)
      P = Point(sphereDeg2cart(lonR[ilon],-lat[ilat],RadEarth))
      (Face_id, iPosFace_id, jPosFace_id) = walk_to_nc(P,start_Face,xw,TransSphereS,RadEarth,Grid)
      start_Face = Face_id
      Inside = InsideFace(P,Grid.Faces[start_Face],Grid)
      if Inside
        iG = Glob[iPosFace_id,jPosFace_id,Face_id]
        Height[iG] += zLevelR[ilon,ilat]
        NumHeight[iG] += 1
      end  
    end
  end
  ExchangeData!(Height,Global.Exchange)
  ExchangeData!(NumHeight,Global.Exchange)
  @. Height /= (NumHeight + 1.e-14)
  @inbounds for iF = 1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = Glob[iP,jP,iF]
        HeightCG[iP,jP,iF] = Height[ind]
      end
    end
  end
  return HeightCG
end

function BoundingBox(Grid)
  MinLonL = 0.0 
  MaxLonL = -180.0
  MinLonR = 180.0  
  MaxLonR = 0.0
  MinLat = 1.e20
  MaxLat = -1.e20
  for i = 1 : Grid.NumNodes
    P = Grid.Nodes[i].P
    (lon, lat) = cart2sphereDeg(P.x,P.y,P.z)  
    if lon >= 0.0
      MinLonR = min(MinLonR, lon)
      MaxLonR = max(MaxLonR, lon)
    else  
      MinLonL = min(MinLonL, lon)
      MaxLonL = max(MaxLonL, lon)
    end  
    MinLat = min(MinLat, lat)
    MaxLat = max(MaxLat, lat)
  end
  return (MinLonL,MaxLonL,MinLonR,MaxLonR,MinLat,MaxLat)
end
