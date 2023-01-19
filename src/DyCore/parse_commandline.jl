function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--Problem"
    help = "Problem name"
    arg_type = String
    default = "Problem"

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

    "--HorLimit"
    help = "Horizontal limiter"
    arg_type = Bool
    default = false

    "--Upwind"
    help = "Vertical diiferencing"
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
    help = "Coriolis parameeterization"
    arg_type = String
    default = "Sphere"

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

#   Domain decomposition
    "--Decomp"
    help = "Domain decomposition method"
    arg_type = String
    default = "Hilbert"

#   Time integration            
    "--SimDays"
     help = "Number of simulation days"
     arg_type = Float64
     default = 10.0

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
    "--GridType"
    help = "Grid"
    arg_type = String
    default = "CubedSphere"

    "--stretch"
    help = "Grid stretching"
    arg_type = Bool
    default = false

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

    "--OrdPoly"
    help = "Order of the local polynoms"
    arg_type = Int
    default = 4

  end
  return parse_args(s)
end
