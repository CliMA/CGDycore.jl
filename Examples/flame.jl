# mpirun -n 1 julia --project Examples/flame.jl
import CGDycore

ENV["CI_PERF_SKIP_RUN"] = true # we only need haskey(ENV, "CI_PERF_SKIP_RUN") == true

try # capture integrator
    include(joinpath(pkgdir(CGDycore), "Examples", "testNHBaroWaveDrySphere.jl"))
catch err
    if err.error !== :exit_profile
        rethrow(err.error)
    end
end

structs = (; U,dtau,CG,Global,Param)

function do_work!(structs)
  (; U,dtau,CG,Global,Param) = structs
  for i=1:20
      CGDycore.RosenbrockSchur!(
        U,
        dtau,
        CGDycore.FcnNHCurlVecI!,
        CGDycore.JacSchur!,
        CG,
        Global,
        Param
      );
  end
  return nothing
end

# integrator is defined
do_work!(structs) # compile first
import Profile
# # Julia 1.7
# Profile.clear_malloc_data()
# Profile.clear()
# prof = Profile.@profile begin
#     do_work!(structs)
# end

# Julia 1.8
Profile.Allocs.clear()
Profile.clear()
prof = Profile.@profile begin
    do_work!(structs)
end

import ProfileCanvas
include("profile_canvas_patch.jl")

if true # interactive
    ProfileCanvas.view(Profile.fetch())
else # non-interactive
    html_file(joinpath(@__DIR__, "flame.html"))
end

