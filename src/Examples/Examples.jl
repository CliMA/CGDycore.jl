module Examples

import ..Grids
import ..Thermodynamics


function InitialProfile!(Model,Problem,Param,Phys)
  # Initial values
  if Problem == "Galewski"
    Profile = Examples.GalewskiExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "BaroWaveDrySphere" || Problem == "BaroWaveHillDrySphere"
    Profile = Examples.BaroWaveExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "SchaerSphericalSphere"
    Profile = Examples.SchaerSphereExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "HeldSuarezDrySphere"
    Profile, Force = Examples.HeldSuarezDryExample()(Param,Phys)
    Model.InitialProfile = Profile
    Model.Force = Force
  elseif Problem == "HeldSuarezMoistSphere"
    Profile, Force, Eddy = Examples.HeldSuarezMoistExample()(Param,Phys)
    Model.InitialProfile = Profile
    Model.Force = Force
    Model.Eddy = Eddy
  end
end

include("parameters.jl")
include("initial.jl")
include("force.jl")
include("topography.jl")
include("PerturbProfile.jl")

end
