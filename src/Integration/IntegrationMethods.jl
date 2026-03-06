mutable struct RosenbrockMethod{FT<:AbstractFloat} <: IntegrationMethod
  name::String
  nStage::Int
  alpha::Array{FT, 2}
  gamma::Array{FT, 2}
  b::Array{FT, 1}
  a::Array{FT, 2}
  c::Array{FT, 2}
  gammaD::FT
  m::Array{FT, 1}
  JacComp::Bool
end

mutable struct RungeKuttaExMethod{FT<:AbstractFloat} <: IntegrationMethod
  name::String
  nStage::Int
  A::Array{FT, 2}
  b::Array{FT, 1}
  JacComp::Bool
end
