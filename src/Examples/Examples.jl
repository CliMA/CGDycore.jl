module Examples

import ..Grids
import ..Thermodynamics


function InitialProfile!(Model,Problem,Param,Phys)
  # Initial values
  if Problem == "GalewskiSphere"
    Profile = Examples.GalewskiExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "HaurwitzSphere"
    Profile = Examples.HaurwitzExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "LinearBlob"
    Profile = Examples.LinearBlob()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "Advection"
    Profile = Examples.Advec()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "AdvectionSphereSpherical"
    Profile = Examples.AdvectionSphereSpherical()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "BaroWaveDrySphere" || Problem == "BaroWaveHillDrySphere"
    Profile = Examples.BaroWaveExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "SchaerSphericalSphere"
    Profile = Examples.SchaerSphereExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro"
    Profile, Force = Examples.HeldSuarezDryExample()(Param,Phys)
    Model.InitialProfile = Profile
    Model.Force = Force
  elseif Problem == "HeldSuarezMoistSphere" || Problem == "HeldSuarezMoistSphereOro"
    Profile, Force = Examples.HeldSuarezMoistExample()(Param,Phys)
    Model.InitialProfile = Profile
    Model.Force = Force
  elseif Problem == "Stratified" || Problem == "HillAgnesiXCart" || Problem == "HillSchaerCart"
    Profile = Examples.StratifiedExample()(Param,Phys)
    Model.InitialProfile = Profile
    @show "Stratified"
  elseif Problem == "WarmBubble2DXCart"
    Profile = Examples.WarmBubbleCartExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "BryanFritschCart"
    ProfileBF = Models.TestRes(Phys)
    Profile = Examples.BryanFritsch(ProfileBF)(Param,Phys)
    Model.InitialProfile = Profile
  end
end

include("parameters.jl")
include("initial.jl")
include("force.jl")
include("topography.jl")
include("InitProfileBryanFritsch.jl")

end
