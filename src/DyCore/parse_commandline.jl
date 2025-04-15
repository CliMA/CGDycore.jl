function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin

    "--Problem"
    help = "Problem name"
    arg_type = String
    default = "Problem"

    "--Discretization"
    help = "Discretization type"
    arg_type = String
    default = "CG"

    "--FluxDG"
    help = "Average flux for DG"
    arg_type = String
    default = "KennedyGruberGrav"

    "--InterfaceFluxDG"
    help = "Interface flux for DG"
    arg_type = String
    default = "KennedyGruberGrav"

    "--Profile"
    help = "Profile for initial onditions"
    arg_type = Bool
    default = false

    "--NumV"
    help = "Number of variables"
    arg_type = Int
    default = 5

    "--NumTr"
    help = "Number of tracer variables"
    arg_type = Int
    default = 0

    "--NumAux"
    help = "Number of auxilary variables"
    arg_type = Int
    default = 1

    "--ProfRho"
    help = "Initial conditions for density"
    arg_type = String
    default = ""

    "--ProfTheta"
    help = "Initial conditions for potential temperature"
    arg_type = String
    default = ""

    "--PertTh"
    help = "Perturb initial conditions for potential temperature"
    arg_type = Bool
    default = false

    "--ProfVel"
    help = "Initial conditions for velocity"
    arg_type = String
    default = ""

    "--ProfVelGeo"
    help = "Initial conditions for geostrophic velocity"
    arg_type = String
    default = ""

    "--ProfVelW"
    help = "Initial conditions for vertical velocity"
    arg_type = String
    default = ""

    "--ProfTr"
    help = "Initial conditions for tracer"
    arg_type = String
    default = ""

    "--ProfpBGrd"
    help = "Initial conditions for background pressure"
    arg_type = String
    default = ""

    "--ProfRhoBGrd"
    help = "Initial conditions for background density"
    arg_type = String
    default = ""

    "--ProfTest"
    help = "Initial conditions for test operators"
    arg_type = String
    default = ""

    "--RhoTPos"
    help = "Position of total water in the tracer list"
    arg_type = Int
    default = 0

    "--RhoVPos"
    help = "Position of water vapor in the tracer list"
    arg_type = Int
    default = 0

    "--RhoCPos"
    help = "Position of cloud water in the tracer list"
    arg_type = Int
    default = 0

    "--RhoIPos"
    help = "Position of cloud ice in the tracer list"
    arg_type = Int
    default = 0

    "--RhoRPos"
    help = "Position of rain water in the tracer list"
    arg_type = Int
    default = 0

    "--TkePos"
    help = "Position of turbulent kinetic energy in the tracer list"
    arg_type = Int
    default = 0

    "--HorLimit"
    help = "Horizontal limiter"
    arg_type = Bool
    default = false

    "--Upwind"
    help = "Vertical differencing"
    arg_type = Bool
    default = true

    "--Damping"
    help = "Rayleigh damping"
    arg_type = Bool
    default = false

    "--Geos"
    help = "Rayleigh damping with geostrphic wind"
    arg_type = Bool
    default = false

    "--Relax"
    help = "Relaxation parameter [1/s] for Rayleigh damping"
    arg_type = Float64
    default = 0.0

    "--StrideDamp"
    help = "Vertical extent [m]  for Rayleigh damping"
    arg_type = Float64
    default = 0.0

    "--Coriolis"
    help = "Coriolis"
    arg_type = Bool
    default = false

    "--CoriolisType"
    help = "Coriolis parameterization"
    arg_type = String
    default = "Shallow"

    "--Buoyancy"
    help = "Buoyancy "
    arg_type = Bool
    default = true

    "--Turbulence"
    help = "Turbulence "
    arg_type = Bool
    default = false

    "--RefProfile"
    help = "RefProfile"
    arg_type = Bool
    default = false

    "--Equation"
    help = "Equation Type, Shallow, Deep"
    arg_type = String
    default = "Shallow"

    "--EDMF"
    help = "EDMF"
    arg_type = Bool
    default = false

    "--NDEDMF"
    help = "Number of Drafts in EDMF"
    arg_type = Int
    default = 0

    "--State"
    help = "Equation Type, ShallowWater, Dry, Moist"
    arg_type = String
    default = "Dry"

    "--Thermo"
    help = "Thermodynamic variable"
    arg_type = String
    default = "PotTemp"

    "--Microphysics"
    help = "Microphysics"
    arg_type = Bool
    default = false

    "--TypeMicrophysics"
    help = "TypeMicrophysics"
    arg_type = String
    default = ""

    "--RelCloud"
    help = "Relaxation parameter [1/s] for cloud microphysics"
    arg_type = Float64
    default = 0.0

    "--Rain"
    help = "Remove cloud water"
    arg_type = Float64
    default = 0.0

    "--Sedimentation"
    help = "Sedimentation"
    arg_type = Bool
    default = false

    "--Source"
    help = "Source"
    arg_type = Bool
    default = false

    "--Forcing"
    help = "Forcing"
    arg_type = Bool
    default = false

    "--JacVerticalDiffusion"
    help = "Jacobian Vertical diffusion"
    arg_type = Bool
    default = false

    "--JacVerticalAdvection"
    help = "Jacobian Vertical advection"
    arg_type = Bool
    default = false

    "--VerticalDiffusion"
    help = "Vertical diffusion"
    arg_type = Bool
    default = false

    "--VerticalDiffusionMom"
    help = "Vertical diffusion for momentum"
    arg_type = Bool
    default = false

    "--SurfaceFlux"
    help = "Surface flux for scalars"
    arg_type = Bool
    default = false

    "--SurfaceFluxMom"
    help = "Surface flux for momentum"
    arg_type = Bool
    default = false

    "--SurfaceScheme"
    help = "Surface scheme"
    arg_type = String
    default = "MOST"

    "--BoundaryWE"
    help = "Boundary type in west east direction"
    arg_type = String
    default = ""

    "--BoundarySN"
    help = "Boundary type in south north direction"
    arg_type = String
    default = ""

    "--BoundaryBT"
    help = "Boundary type in vertical direction"
    arg_type = String
    default = ""

    "--Curl"
    help = "Form of momentum transport, vector invariant or advective"
    arg_type = Bool
    default = true

    "--ModelType"
    help = "Form of momentum transport, vector invariant or advective or conservative"
    arg_type = String
    default = "VectorInvariant"


#   Domain decomposition
    "--Decomp"
    help = "Domain decomposition method"
    arg_type = String
    default = "EqualArea"

#   Time integration            
    "--SimDays"
     help = "Number of simulation days"
     arg_type = Int
     default = 0

    "--SimHours"
     help = "Number of simulation hours"
     arg_type = Int
     default = 0

    "--SimMinutes"
     help = "Number of simulation minutes"
     arg_type = Int
     default = 0

    "--SimSeconds"
     help = "Number of simulation seconds"
     arg_type = Int
     default = 0

    "--SimTime"
     help = "Simulation time"
     arg_type = Float64
     default = 0.0

    "--IntMethod"
    help = "Integration  method"
    arg_type = String
    default = "Rosenbrock"

    "--Table"
    help = "Butcher tableaux"
    arg_type = String
    default = "SSP-Knoth"

    "--dtau"
    help = "Time step"
    arg_type = Float64
    default = 100.0

#   Orography
    "--TopoS"
    help = "Orography type"
    arg_type = String
    default = ""

#   Grid
    "--GridForm"
    help = "Cartesian or Spherical"
    arg_type = String
    default = "Spherical"

    "--RadEarth"
    help = "Radius of sphere"
    arg_type = Float64
    default = 0.0 

    "--ScaleFactor"
    help = "ScaleFactor for EarthRadius"
    arg_type = Float64
    default = 0.0 

    "--GridType"
    help = "Grid"
    arg_type = String
    default = "CubedSphere"

    "--Stretch"
    help = "Grid stretching"
    arg_type = Bool
    default = false

    "--StretchType"
    help = "Type of grid stretching, ICON, Exp "
    arg_type = String
    default = ""

    "--AdaptGridType"
    help = "Type of boundary following coordinate"
    arg_type = String
    default = "Sleve"

    "--nx"
    help = "Number of horizontal grid cell in x-direction"
    arg_type = Int
    default = 1

    "--ny"
    help = "Number of horizontal grid cell in y-direction"
    arg_type = Int
    default = 1

    "--nz"
    help = "Number of vertical levels"
    arg_type = Int
    default = 1

    "--H"
    help = "Domain height in [m]"
    arg_type = Float64
    default = 1000.0

    "--nPanel"
    help = "Number of grid cells per panel for a cubed sphere grid"
    arg_type = Int
    default = 4

    "--ns"
    help = "Number of grid cells per healpix"
    arg_type = Int
    default = 60

    "--RefineLevel"
    help = "Number of refinement levels for a triangular grid"
    arg_type = Int
    default = 0

    "--nLon"
    help = "Number of grid cells in longitudidal direction"
    arg_type = Int
    default = 1

    "--nLat"
    help = "Number of grid cells in latitudidal direction"
    arg_type = Int
    default = 1

    "--LatB"
    help = "Pol cap"
    arg_type = Float64
    default = 8 / 9 * 0.5 * pi

    "--Lx"
    help = "Length in [m] in x-direction"
    arg_type = Float64
    default = 1000.0

    "--Ly"
    help = "Length in [m] in y-direction"
    arg_type = Float64
    default = 1000.0

    "--x0"
    help = "Initial point in [m] in x-direction"
    arg_type = Float64
    default = 1000.0

    "--y0"
    help = "Initial point in [m] in y-direction"
    arg_type = Float64
    default = 1000.0

    "--OrdPoly"
    help = "Order of the local polynoms"
    arg_type = Int
    default = 4

    "--OrdPolyZ"
    help = "Order of the local polynoms in vertical direction"
    arg_type = Int
    default = 0

    "--HyperVisc"
    help = "HyperViscosity"
    arg_type = Bool
    default = false

    "--HyperDCurl"
    help = "HyperDCurl"
    arg_type = Float64
    default = 0.0

    "--HyperDGrad"
    help = "HyperDGrad"
    arg_type = Float64
    default = 0.0

    "--HyperDRhoDiv"
    help = "HyperDRhoDiv"
    arg_type = Float64
    default = 0.0

    "--HyperDDiv"
    help = "HyperDDiv"
    arg_type = Float64
    default = 0.0

    "--HyperDDivW"
    help = "HyperDDivW"
    arg_type = Float64
    default = 0.0

    "--P1"
    help = "P1"
    arg_type = Float64
    default = 0.0

    "--P2"
    help = "P2"
    arg_type = Float64
    default = 0.0

    "--P3"
    help = "P3"
    arg_type = Float64
    default = 0.0

    "--P4"
    help = "P4"
    arg_type = Float64
    default = 0.0

#   Output    

    "--vtkFileName"
     help = "File mame of vtk output"
     arg_type = String
     default = "Output"

    "--OrdPrint"
     help = "Subcell printing"
     arg_type = Int
     default = 0

    "--PrintDays"
     help = "Number of print days"
     arg_type = Int
     default = 0

    "--PrintHours"
     help = "Number of print hours"
     arg_type = Int
     default = 0

    "--PrintMinutes"
     help = "Number of print minutes"
     arg_type = Int
     default = 0

    "--PrintSeconds"
     help = "Number of print seconds"
     arg_type = Int
     default = 0

    "--PrintTime"
     help = "Print interval"
     arg_type = Float64
     default = 0.0

    "--PrintStartTime"
     help = "Print start"
     arg_type = Float64
     default = 0.0

    "--StartAverageDays"
     help = "Start statistics"
     arg_type = Int
     default = -1

    "--Flat"
     help = "Output as sphere or cube"
     arg_type = Bool
     default = false

     "--RefineOutput"
     help = "Refining Output Finite Elements"
     arg_type = Int
     default = 0

    "--Device"
     help = "CPU or GPU"
     arg_type = String
     default = "CPU"

    "--GPUType"
     help = "Type of GPU, CUDA, AMD, Metal"
     arg_type = String
     default = "CUDA"

    "--FloatTypeBackend"
     help = "Type of Float for Backend"
     arg_type = String
     default = "Float32"

    "--NumberThreadGPU"
     help = "Number of threads for GPU"
     arg_type = Int
     default = 256

    "--NumberThreadTriGPU"
     help = "Number of threads for linear vertical solver in GPU mode"
     arg_type = Int
     default = 32

  # Finite elements     

    "--OrderFEM"
     help = "Order of finite elements"
     arg_type = Int
     default = 0

  end
  return parse_args(s)
end
