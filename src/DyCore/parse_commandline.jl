function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin

    "--Problem"
    help = "Problem name"
    arg_type = String
    default = "Problem"

    "--NumV"
    help = "Number of variables"
    arg_type = Int
    default = 5

    "--NumTr"
    help = "Number of tracer variables"
    arg_type = Int
    default = 0

    "--ProfRho"
    help = "Initial conditions for density"
    arg_type = String
    default = ""

    "--ProfTheta"
    help = "Initial conditions for potential temperature"
    arg_type = String
    default = ""

    "--ProfVel"
    help = "Initial conditions for velocity"
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
    default = "Sphere"

    "--Buoyancy"
    help = "Buoyancy "
    arg_type = Bool
    default = true

    "--RefProfile"
    help = "RefProfile"
    arg_type = Bool
    default = false

    "--Equation"
    help = "Equation Type, Shallow, Compressible, CompressibleMoist"
    arg_type = String
    default = "Compressible"

    "--Microphysics"
    help = "Microphysics"
    arg_type = Bool
    default = false

    "--Source"
    help = "Source"
    arg_type = Bool
    default = false

    "--VerticalDiffusion"
    help = "Vertical diffusion"
    arg_type = Bool
    default = false

    "--SurfaceFlux"
    help = "Surface flux"
    arg_type = Bool
    default = false

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

    "--Thermo"
    help = "Type of thermodynmic variable, total energy, internal energy, potential temperature"
    arg_type = String
    default = ""

    "--Curl"
    help = "Form of momentum transport, vector invariant or advective"
    arg_type = Bool
    default = true


#   Domain decomposition
    "--Decomp"
    help = "Domain decomposition method"
    arg_type = String
    default = "Hilbert"

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
    "--RadEarth"
    help = "Radius of sphere"
    arg_type = Float64
    default = 0.0 

    "--GridType"
    help = "Grid"
    arg_type = String
    default = "CubedSphere"

    "--stretch"
    help = "Grid stretching"
    arg_type = Bool
    default = false

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

    "--HyperDDiv"
    help = "HyperDDiv"
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

#   Output    

    "--vtkFileName"
     help = "File mame of vtk output"
     arg_type = String
     default = "Output"

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

    "--Flat"
     help = "Output as sphere or cube"
     arg_type = Bool
     default = true

  end
  return parse_args(s)
end
