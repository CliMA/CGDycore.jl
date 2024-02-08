abstract type QuadType end


mutable struct QuadRule{FT<:AbstractFloat,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray} <: QuadType
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

function QuadRule{FT}(type::Grids.Tri,backend,n) where FT<:AbstractFloat

  if n == 1
    NumQuad = 1  
    WeightsCPU = zeros(Float64,NumQuad)
    PointsCPU = zeros(Float64,NumQuad,2)
    Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
    Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
    WeightsCPU[1] = 2
    PointsCPU[1,1] = -1/3
    PointsCPU[1,2] = -1/3
  elseif n == 2
    NumQuad = 3  
    WeightsCPU = zeros(Float64,NumQuad)
    PointsCPU = zeros(Float64,NumQuad,2)
    Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
    Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
    WeightsCPU[1] = 2/3
    WeightsCPU[2] = 2/3
    WeightsCPU[3] = 2/3
    PointsCPU[1,1] = -2/3
    PointsCPU[1,2] = 1/3
    PointsCPU[2,1] = 1/3
    PointsCPU[2,2] = -2/3
    PointsCPU[3,1] = -2/3
    PointsCPU[3,2] = -2/3
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
