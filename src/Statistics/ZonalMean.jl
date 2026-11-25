# Zonal and time means on latitude bands (--ZonalMeanHours), sampled on the host from
# StartAverageDays on. Level fields are volume weighted with the assembled nodal mass, surface
# fields area weighted with the mass of the lowest layer divided by its thickness.
# Precipitation is the condensation rate of the large-scale condensation with immediate
# fall-out (Rain = 1): RelCloud * max(rho_v - rho_vs(T), 0), summed over the column, or with
# --TypeMicrophysics="Adjustment" the condensed mass of the last step (CGSEM.Condensation).
const ZonalLevelNames = ["u","v","w","T","theta","qv","p","uu","vv","uv","vT"]
const ZonalSurfNames = ["precip","evap","shf","ustar","TS","cwv"]

mutable struct ZonalMeanStruct
  nLat::Int
  Interval::Float64
  StartTime::Float64
  NextTime::Float64
  Bin::Vector{Int}
  W::Matrix{Float64}
  WS::Vector{Float64}
  Sum::Array{Float64,3}
  SumW::Matrix{Float64}
  SumS::Matrix{Float64}
  SumWS::Vector{Float64}
  nSample::Int
  FileName::String
end

function ZonalMeanStruct(Metric,nz,Hours,StartDays,FileName;nLat=90)
  xS = Array(Metric.xS)
  MH = Array(Metric.M)
  dzH = Array(Metric.dz)
  NumG = size(xS,2)
  Bin = [clamp(floor(Int,(xS[2,ic] + pi / 2) / pi * nLat) + 1,1,nLat) for ic = 1 : NumG]
  W = zeros(nz,NumG)
  WS = zeros(NumG)
  for ic = 1 : NumG
    for iz = 1 : nz
      W[iz,ic] = MH[iz,ic,1] + MH[iz,ic,2]
    end
    WS[ic] = W[1,ic] / dzH[1,ic]
  end
  StartTime = 86400.0 * max(StartDays,0)
  return ZonalMeanStruct(nLat,3600.0 * Hours,StartTime,StartTime,Bin,W,WS,
    zeros(nz,nLat,length(ZonalLevelNames)),zeros(nz,nLat),
    zeros(nLat,length(ZonalSurfNames)),zeros(nLat),0,FileName)
end

function ZonalMeanSample!(ZM::ZonalMeanStruct,time,U,p,SurfaceData,Metric,Phys,Model)
  if time + 1.0e-6 < ZM.NextTime
    return
  end
  ZM.NextTime += ZM.Interval
  UH = Array(U)
  pH = Array(p)
  dzH = Array(Metric.dz)
  SD = Array(SurfaceData)
  PrecipH = Model.Condensation !== nothing ? Array(Model.Condensation.Precip) : nothing
  RhoVPos = Model.RhoVPos
  RhoCPos = Model.RhoCPos
  RelCloud = Model.RelCloud
  nz = size(UH,2)
  kappa = Phys.Rd / Phys.Cpd
  for ic in axes(UH,3)
    ib = ZM.Bin[ic]
    Precip = 0.0
    CWV = 0.0
    for iz = 1 : nz
      Rho = UH[1,iz,ic,1]
      RhoV = RhoVPos > 0 ? UH[1,iz,ic,RhoVPos] : 0.0
      RhoC = RhoCPos > 0 ? UH[1,iz,ic,RhoCPos] : 0.0
      RhoD = Rho - RhoV - RhoC
      pL = pH[1,iz,ic]
      T = pL / (Phys.Rd * RhoD + Phys.Rv * RhoV)
      u = UH[1,iz,ic,2]
      v = UH[1,iz,ic,3]
      w = UH[1,iz,ic,4]
      Val = (u,v,w,T,UH[1,iz,ic,5] / Rho,RhoV / Rho,pL,u * u,v * v,u * v,v * T)
      Wi = ZM.W[iz,ic]
      for iv in eachindex(Val)
        ZM.Sum[iz,ib,iv] += Wi * Val[iv]
      end
      ZM.SumW[iz,ib] += Wi
      if RhoVPos > 0
        if PrecipH === nothing
          RhoVS = Thermodynamics.fpws(T,Phys) / (Phys.Rv * T)
          Precip += RelCloud * max(RhoV - RhoVS,0.0) * dzH[iz,ic]
        end
        CWV += RhoV * dzH[iz,ic]
      end
    end
    if PrecipH !== nothing
      Precip = PrecipH[ic]
    end
    # surface fluxes as in Surfaces.SurfaceFlux
    Rho1 = UH[1,1,ic,1]
    RhoV1 = RhoVPos > 0 ? UH[1,1,ic,RhoVPos] : 0.0
    T1 = pH[1,1,ic] / (Phys.Rd * (Rho1 - RhoV1) + Phys.Rv * RhoV1)
    pS = pH[1,1,ic] + Rho1 * Phys.Grav * 0.5 * dzH[1,ic]
    TS = SD[Surfaces.TSurfPos,ic]
    Ua = SD[Surfaces.UaPos,ic]
    Evap = RhoVPos > 0 ? SD[Surfaces.CHPos,ic] * Ua * (SD[Surfaces.RhoVSurfPos,ic] - RhoV1) : 0.0
    SHF = Phys.Cpd * Rho1 * SD[Surfaces.CTPos,ic] * Ua * (TS * (pH[1,1,ic] / pS)^kappa - T1)
    ValS = (Precip,Evap,SHF,SD[Surfaces.uStarPos,ic],TS,CWV)
    for iv in eachindex(ValS)
      ZM.SumS[ib,iv] += ZM.WS[ic] * ValS[iv]
    end
    ZM.SumWS[ib] += ZM.WS[ic]
  end
  ZM.nSample += 1
end

# binary Float64: level fields (nz,nLat,nVar), then surface fields (nLat,nSurf); header in .txt
function ZonalMeanWrite(ZM::ZonalMeanStruct,time)
  if ZM.nSample == 0
    return
  end
  nz, nLat, nVar = size(ZM.Sum)
  Mean = zeros(nz,nLat,nVar)
  for iv = 1 : nVar, ib = 1 : nLat, iz = 1 : nz
    Mean[iz,ib,iv] = ZM.Sum[iz,ib,iv] / max(ZM.SumW[iz,ib],eps())
  end
  MeanS = zeros(nLat,length(ZonalSurfNames))
  for iv in axes(MeanS,2), ib = 1 : nLat
    MeanS[ib,iv] = ZM.SumS[ib,iv] / max(ZM.SumWS[ib],eps())
  end
  open(ZM.FileName * ".bin","w") do io
    write(io,Mean)
    write(io,MeanS)
  end
  open(ZM.FileName * ".txt","w") do io
    println(io,"nz $nz nLat $nLat nVar $nVar nSurf $(length(ZonalSurfNames)) samples $(ZM.nSample) ",
      "start $(ZM.StartTime) end $time")
    println(io,"level ",join(ZonalLevelNames," "))
    println(io,"surface ",join(ZonalSurfNames," "))
    println(io,"lat bins: equal width from -90 to 90 degrees, level fields at model levels (not pressure)")
  end
end

# extrema of the state with position (lat, lon in degrees, level), printed with the outputs
function ExtremaPrint(time,U,p,Metric,Phys,Model)
  UH = Array(U)
  pH = Array(p)
  xS = Array(Metric.xS)
  RhoVPos = Model.RhoVPos
  nz = size(UH,2)
  Ext = Dict{String,Tuple{Float64,Int,Int}}()
  function upd!(name,val,iz,ic,larger)
    if !haskey(Ext,name) || (larger ? val > Ext[name][1] : val < Ext[name][1]) || isnan(val)
      Ext[name] = (val,iz,ic)
    end
  end
  for ic in axes(UH,3), iz = 1 : nz
    Rho = UH[1,iz,ic,1]
    RhoV = RhoVPos > 0 ? UH[1,iz,ic,RhoVPos] : 0.0
    T = pH[1,iz,ic] / (Phys.Rd * (Rho - RhoV) + Phys.Rv * RhoV)
    upd!("max|w|",abs(UH[1,iz,ic,4]),iz,ic,true)
    upd!("max|u|",hypot(UH[1,iz,ic,2],UH[1,iz,ic,3]),iz,ic,true)
    upd!("minT",T,iz,ic,false)
    upd!("maxT",T,iz,ic,true)
    if RhoVPos > 0
      upd!("maxqv",RhoV / Rho,iz,ic,true)
      upd!("maxRH",RhoV * Phys.Rv * T / Thermodynamics.fpws(T,Phys),iz,ic,true)
    end
  end
  s = "Extrema t=$(round(time / 86400,digits=3)) d:"
  for name in sort(collect(keys(Ext)))
    val, iz, ic = Ext[name]
    s *= " $name $(round(val,sigdigits=4)) (lat $(round(rad2deg(xS[2,ic]),digits=1)) " *
      "lon $(round(rad2deg(xS[1,ic]),digits=1)) iz $iz)"
  end
  println(s)
  flush(stdout)
end
