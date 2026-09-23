module Parallels

using MPI
using StaticArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll

include("Exchange.jl")
include("Hilbert.jl")
include("EqualAreaPartitioner.jl")

end
