abstract type QuadType end


mutable struct QuadRule{FT<:AbstractFloat,
                        IT1<:AbstractArray,
                        IT2<:AbstractArray} <: QuadType
    NumQuad::Int
    Weights::IT1
    Points::IT2 
end

function QuadRule{FT}(type::Grids.QuadPrimal,backend,n) where FT<:AbstractFloat
  if n == 1
    NumQuad = 4
    x, w = gaussradau(2)
    WeightsCPU = zeros(NumQuad)
    PointsCPU = zeros(NumQuad,2)
    Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
    Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
    ii = 1
    for j = 1 : 2
      for i = 1 : 2
        WeightsCPU[ii] = w[i] * w[j]
        PointsCPU[ii,1] = x[i]
        PointsCPU[ii,2] = x[j]
        ii += 1
      end
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


function QuadRule{FT}(type::Grids.QuadDual,backend,n) where FT<:AbstractFloat
  if n == 1
    NumQuad = 4  
    x, w = gaussradau(2)
    x .= -x
    WeightsCPU = zeros(NumQuad)
    PointsCPU = zeros(NumQuad,2)
    Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
    Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
    ii = 1
    for j = 1 : 2
      for i = 1 : 2
        WeightsCPU[ii] = w[i] * w[j]
        PointsCPU[ii,1] = x[i]
        PointsCPU[ii,2] = x[j]
        ii += 1
      end
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

function QuadRule{FT}(type::Grids.Line,backend,n) where FT<:AbstractFloat
  x, w = gausslobatto(n+1)
  NumQuad = (n+1)
  WeightsCPU = zeros(Float64,NumQuad)
  PointsCPU = zeros(Float64,NumQuad)
  Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
  Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
  ii = 1
  for i = 1 : n+1
    WeightsCPU[ii] = w[i]
    PointsCPU[ii,1]=x[i]
    ii += 1
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


function QuadRule{FT}(type::Grids.Quad,backend,n) where FT<:AbstractFloat
    x, w = gausslobatto(n+1)
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
  elseif n == 3
    NumQuad = 4 
    WeightsCPU = zeros(Float64,NumQuad)
    PointsCPU = zeros(Float64,NumQuad,2)
    Weights = KernelAbstractions.zeros(backend,FT,NumQuad)
    Points = KernelAbstractions.zeros(backend,FT,NumQuad,2)
    WeightsCPU[1] = -27/96
    WeightsCPU[2] = 25/96
    WeightsCPU[3] = 25/96
    WeightsCPU[4] = 25/96
    PointsCPU[1,1] = 1/3
    PointsCPU[1,2] = 1/3
    PointsCPU[2,1] = -3/5
    PointsCPU[2,2] = -3/5
    PointsCPU[3,1] = 3/5
    PointsCPU[3,2] = -3/5
    PointsCPU[4,1] = -3/5
    PointsCPU[4,2] = 3/5
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
