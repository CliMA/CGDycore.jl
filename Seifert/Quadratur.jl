abstract type QuadType end


mutable struct QuadRule{FT<:AbstractFloat,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray} <: FiniteElement
NumQuad::Int
Weights::IT1
Points::IT2 
end

function QuadRule{FT}(type::Grids.Quad,backend,n) where FT<:AbstractFloat
    w,x = DG.GaussLobattoQuad(n)
    NumQuad = (n+1)*(n+1)
    WeightsCPU = zeros(Float64,NumQuad)
    PointsCPU = zeros(Float64,NumQuad,2)
    Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
    Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
    ii = 1
    for j = 1 : n+1
        for i = 1 : n+1
            WeightsCPU[ii] = w[i]*w[j]
            PointsCPU[ii,1]=x[i]
            PointsCPU[ii,2]=x[j] 
            ii += 1
        end    
    end
    copyto!(Weights,WeightsCPU)
    copyto!(Points,PointsCPU)
    return QuadRule{FT,
                    typeof(Weights),
                    typeof(Points)}( 
    NumQuad,
    Weights,
    Points,
      )
end
