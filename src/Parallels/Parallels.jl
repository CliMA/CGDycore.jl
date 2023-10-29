module Parallels

using MPI
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace

include("Exchange.jl")
include("Hilbert.jl")
include("EqualAreaPartitioner.jl")

end
