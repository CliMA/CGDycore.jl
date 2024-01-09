module FiniteElements

using Polynomials
using SpecialPolynomials

abstract type FiniteElement end

struct NodalElement <: FiniteElement 
  NumBases::Int
  NodalBases::Array{Polynomial,1}
  NodalPoints::Array{Float64,2}
  Diff1Matrix::Array{Float64,2}
end

struct W1 <: FiniteElement
  Order::Int
  phiv::Array{Polynomial,3}
  phie::Array{Polynomial,3}
end

function W1(Order)
  phiv = Array{Polynomial,3}(undef,1,2,3)
  phiv[1,1,1] = Polynomial([0.5,-0.5])
  phiv[1,1,2] = Polynomial([0.5,0.5])
  phiv[1,1,3] = Polynomial([1])
  phiv[1,2,1] = Polynomial([0.5,-0.5])
  phiv[1,2,2] = Polynomial([0.5,-0.5])
  phiv[1,2,3] = Polynomial([0.5,0.5])
  if Order > 3
    phie = Array{Polynomial,3}(undef,Order-2,2,3)
    for m = 1 : Order - 2
      phie[m,1,1] = Polynomial([0.5,0.5]) * 
                    Polynomial([0.5,-0.5])
    end  
  else
    phie = Array{Polynomial,3}(undef,0,2,3)
  end  
      
  return W1(
    Order,
    phiv,
    phie,
  )  
end

end
