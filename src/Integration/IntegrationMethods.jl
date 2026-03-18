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

mutable struct LSRungeKuttaMethod{FT<:AbstractFloat} <: IntegrationMethod
  name::String
  nStage::Int
  a::Array{FT, 1}
  b::Array{FT, 1}
  c::Array{FT, 1}
  JacComp::Bool
end



mutable struct MISMethod{FT<:AbstractFloat} <: IntegrationMethod
  name::String
  nStage::Int
  beta::Array{FT, 2}
  alpha::Array{FT, 2}
  gamma::Array{FT, 2}
  d::Array{FT, 1}
  FastMethod::IntegrationMethod
  JacComp::Bool
end

mutable struct CacheMISStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  dZn::AT4
  fV::AT4
  Sdu::AT5
  Yn::AT5
  CacheFast
end

function Cache(backend,FT,IntMethod::MISMethod,FE,M,nz,NumV)
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)

  dZn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  Sdu = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  Yn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage+1)

  CacheFast = Cache(backend,FT,IntMethod.FastMethod,FE,M,nz,NumV)

  return CacheMISStruct{FT,
                     typeof(Vn),
                     typeof(Sdu)}(
    Vn,
    dZn,
    fV,
    Sdu,
    Yn,
    CacheFast,
  )
end


