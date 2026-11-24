# Boundary-layer test runs of CGSEM without VTK output.
# usage: julia --project Examples/PBLRun.jl <DriverCG arguments>
# environment:
#   PBLTAG   name of the run (default "pbl")
#   PBLDIR   output directory (default "pbl_out")
#   PBLOUT   output interval [s] (default 600)
#   PBLMODE  "column": profiles of the first column as text,
#            "slice":  all columns in the plane y = y0 as Float32 binary
# Column file columns: z u v w theta tke K (w and K at the upper interface of the layer).
# Slice files: <tag>_<time>.bin holds Float32 [var, iz, column] with var = x z u w theta tke K.
src = read("Examples/DriverCG.jl", String)
include_string(Main, src[1:findfirst("Integration.TimeStepper(", src)[1]-1], "DriverCG.jl")

FT = eltype(U)
Tag = get(ENV, "PBLTAG", "pbl")
Dir = get(ENV, "PBLDIR", "pbl_out")
OutInt = parse(Float64, get(ENV, "PBLOUT", "600"))
Mode = get(ENV, "PBLMODE", "column")
mkpath(Dir)
M = CG.Mz
nz = Grid.nz
CacheInt = Integration.Cache(backend, FT, MethodInt, CG, M, nz, size(U, 4))
CacheAux = Integration.CacheAuxStruct(backend, FT, CG, M, nz, Global.Model, Grid)
JCache = MethodInt.JacComp ? Integration.CacheJac(backend, FT, M, nz, Global.Model, CG) : nothing
Global.Model.GPAuxPos > 0 && DGSEM.GeoPot(CacheAux.Aux, CG, Metric, Exchange, Global)
SimTimeAll = 86400 * SimDays + 3600 * SimHours + 60 * SimMinutes + SimSeconds + SimTime
nStep = round(Int, SimTimeAll / dtau)
nOut = max(1, round(Int, OutInt / dtau))
TkePosR = Global.Model.TkePos

zPH = Array(Metric.zP)
GlobH = Array(CG.Glob)
XH = Array(Metric.X)
# columns of the slice: lowest y, sorted by x
SliceCols = Int[]
SliceX = Float64[]
let y0 = minimum(XH[:, 1, 2, 1, :]), seen = Set{Int}()
  for IF in axes(GlobH, 2), ID in axes(GlobH, 1)
    ind = GlobH[ID, IF]
    if abs(XH[ID, 1, 2, 1, IF] - y0) < 1.0e-6 * max(1.0, abs(y0)) + 1.0e-3 && !(ind in seen)
      push!(seen, ind); push!(SliceCols, ind); push!(SliceX, XH[ID, 1, 1, 1, IF])
    end
  end
end
perm = sortperm(SliceX)
SliceCols = SliceCols[perm]; SliceX = SliceX[perm]
if Model.LateralBC !== nothing
  # without the ghost columns east of the outflow
  keep = SliceX .<= Model.LateralBC.xOut + 1.0e-3
  SliceCols = SliceCols[keep]; SliceX = SliceX[keep]
end

function output(t)
  Uh = Array(U)
  KVh = Array(CacheAux.KV)
  Th(iz, ic) = Uh[1, iz, ic, 5] / Uh[1, iz, ic, 1]
  Tke(iz, ic) = TkePosR > 0 ? Uh[1, iz, ic, TkePosR] / Uh[1, iz, ic, 1] : 0.0
  Kin(iz, ic) = KVh[1, iz, ic] / (iz < nz ? 0.5 * (Uh[1, iz, ic, 1] + Uh[1, iz+1, ic, 1]) : Uh[1, iz, ic, 1])
  ts = lpad(round(Int, t), 6, "0")
  if Mode == "column"
    ic = 1
    open(joinpath(Dir, "$(Tag)_$(ts).txt"), "w") do io
      zs = zPH[1, ic] - 0.5 * Array(Metric.dz)[1, ic]
      println(io, "# t = $t s, z above ground; w and K at the upper interface")
      if Global.Model.SurfaceFlux || Global.Model.SurfaceFluxMom
        SD = Array(Global.SurfaceData.Data)
        println(io, "# uStar ", SD[Surfaces.uStarPos, ic], " Ua ", SD[Surfaces.UaPos, ic],
          " CM ", SD[Surfaces.CMPos, ic], " CT ", SD[Surfaces.CTPos, ic])
      end
      for iz in 1:nz
        println(io, zPH[iz, ic] - zs, " ", Uh[1, iz, ic, 2], " ", Uh[1, iz, ic, 3], " ", Uh[1, iz, ic, 4],
          " ", Th(iz, ic), " ", Tke(iz, ic), " ", Kin(iz, ic))
      end
    end
  else
    A = zeros(Float32, 7, nz, length(SliceCols))
    for (j, ic) in enumerate(SliceCols), iz in 1:nz
      A[:, iz, j] .= (SliceX[j], zPH[iz, ic], Uh[1, iz, ic, 2], Uh[1, iz, ic, 4], Th(iz, ic), Tke(iz, ic), Kin(iz, ic))
    end
    write(joinpath(Dir, "$(Tag)_$(ts).bin"), A)
    open(joinpath(Dir, "$(Tag).txt"), "w") do io
      println(io, "nvar 7 nz $nz ncol $(length(SliceCols)) vars x z u w theta tke K")
    end
  end
end

step!() = (Integration.TimeIntegration!(MethodInt, U, dt, Fcn, CacheAux, CGSEM.JacGPU!, CG, Metric, Phys,
  CacheInt, JCache, Exchange, Global, Param, Model.Equation); Integration.ProjectState!(U, CG, Metric, Global))

output(0.0)
tWall = @elapsed for i in 1:nStep
  step!()
  if mod(i, nOut) == 0 || i == nStep
    KernelAbstractions.synchronize(backend)
    finite = all(isfinite, U)
    println("PBL $Tag step $i t $(i * dtau) s finite $finite")
    output(i * dtau)
    finite || break
  end
end
println("PBL $Tag $(eltype(U)) steps $nStep dtau $dtau wall $(round(tWall, digits=1)) s")
