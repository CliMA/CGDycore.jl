module FiniteElements

using Polynomials

abstract type FiniteElement end

struct NodalElement <: FiniteElement 
  NumBases::Int
  NodalBases::Array{Polynomial,1}
  NodalPoints::Array{Float64,2}
  Diff1Matrix::Array{Float64,2}
end


end
