module Sources

import ..Parameters as P
import ..Examples
import ..Grids

using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras

include("damping.jl")
include("coriolis.jl")
include("gravitation.jl")
include("kernel.jl")

end
