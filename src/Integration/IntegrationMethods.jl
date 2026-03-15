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

mutable struct IMEXDirkMethod{FT<:AbstractFloat} <: IntegrationMethod
  name::String
  nStage::Int
  AE::Array{FT, 2}
  bE::Array{FT, 1}
  AI::Array{FT, 2}
  bI::Array{FT, 1}
  D::Array{FT, 2}
  E::Array{FT, 2}
  d::Array{FT, 1}
  e::Array{FT, 1}
  JacComp::Bool
end

mutable struct LinIMEXMethod{FT<:AbstractFloat} <: IntegrationMethod
  name::String
  nStage::Int
  AE::Array{FT, 2}
  bE::Array{FT, 1}
  AI::Array{FT, 2}
  bI::Array{FT, 1}
  D::Array{FT, 2}
  E::Array{FT, 2}
  d::Array{FT, 1}
  e::Array{FT, 1}
# AHat::Array{Float64, 2}
# BHat::Array{Float64, 1}
# A::Array{Float64, 2}
# B::Array{Float64, 1}
  JacComp::Bool
end
