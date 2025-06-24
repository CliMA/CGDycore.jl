module Examples

import ..Grids
import ..Thermodynamics

using NLsolve
using KernelAbstractions

include("parameters.jl")
include("initial.jl")
include("force.jl")
include("topography.jl")
include("InitProfileBryanFritsch.jl")

function InitialProfile!(backend,FTB,Model,Problem,Param,Phys)
  # Initial values
  if Problem == "GalewskySphere"
    Profile = Examples.GalewskyExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "HaurwitzSphere"
    Profile = Examples.HaurwitzExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "ModonCollision"
    Profile = Examples.ModonCollisionExample()(Param,Phys)
    Model.InitialProfile = Profile  
  elseif Problem == "BickleyJet"
    Profile = Examples.BickleyJetExample()(Param,Phys)
    Model.InitialProfile = Profile  
  elseif Problem == "LinearGravity"
    Profile = Examples.LinearGravityExample()(Param,Phys)
    Model.InitialProfile = Profile  
  elseif Problem == "InertiaGravityShortCart" || Problem == "InertiaGravityLongCart"
    Profile = Examples.InertiaGravityExample()(Param,Phys)
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
    Profile = Examples.BaroWaveDryExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == Problem == "BaroWaveMoistSphere"
    Profile = Examples.BaroWaveMoistExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "SchaerSphericalSphere"
    Profile = Examples.SchaerSphereExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "GapSphere"
    Profile = Examples.GapSphereExample()(Param,Phys)
    Model.InitialProfile = Profile
  elseif Problem == "HeldSuarezDrySphere" || Problem == "HeldSuarezDrySphereOro" ||
    Problem == "FriersonSphere"
    Profile, Force = Examples.HeldSuarezDryExample()(Param,Phys)
    Model.InitialProfile = Profile
    Model.Force = Force
  elseif Problem == "HeldSuarezMoistSphere" || Problem == "HeldSuarezMoistSphereOro"
    _, Force = Examples.HeldSuarezMoistExample()(Param,Phys)
#   Model.InitialProfile = Profile
    Profile = Examples.StratifiedSphereExample()(Param,Phys)
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
    ProfileBFCPU = TestRes(Phys)
    ProfileBF = KernelAbstractions.zeros(backend,FTB,size(ProfileBFCPU))
    copyto!(ProfileBF,ProfileBFCPU)
    Profile = Examples.BryanFritsch()(Param,Phys,ProfileBF)
    Model.InitialProfile = Profile
  end
end


end
