function QuadRule(type::Grids.QuadPrimal,n)
  if n == 1
    NumQuad = 4
    x, w = gaussradau(2)
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    ii = 1
    for j = 1 : 2
      for i = 1 : 2
        Weights[ii] = w[i] * w[j]
        Points[ii,1] = x[i]
        Points[ii,2] = x[j]
        ii += 1
      end
    end
  end

  return NumQuad, Weights, Points
end


function QuadRule(type::Grids.QuadDual,n)
  if n == 1
    NumQuad = 4  
    x, w = gaussradau(2)
    x .= -x
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    ii = 1
    for j = 1 : 2
      for i = 1 : 2
        Weights[ii] = w[i] * w[j]
        Points[ii,1] = x[i]
        Points[ii,2] = x[j]
        ii += 1
      end
    end
  end  

  return NumQuad, Weights, Points
end

function QuadRule(type::Grids.Line,n)
  x, w = gausslobatto(n+1)
  NumQuad = (n+1)
  Weights = zeros(NumQuad)
  Points = zeros(NumQuad)
  ii = 1
  for i = 1 : n+1
    Weights[ii] = w[i]
    Points[ii,1]=x[i]
    ii += 1
  end
  return NumQuad, Weights, Points
end


function QuadRule(type::Grids.Quad,n)
    x, w = gausslobatto(n+1)
    NumQuad = (n+1)*(n+1)
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    ii = 1
    for j = 1 : n + 1
        for i = 1 : n + 1
            Weights[ii] = w[i]*w[j]
            Points[ii,1]=x[i]
            Points[ii,2]=x[j] 
            ii += 1
        end    
    end
  return NumQuad, Weights, Points
end

function QuadRule(type::Grids.Tri,n)

  if n == 1
    NumQuad = 1  
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = 2
    Points[1,1] = -1/3
    Points[1,2] = -1/3
  elseif n == 2
    NumQuad = 3  
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = 2/3
    Weights[2] = 2/3
    Weights[3] = 2/3
    Points[1,1] = -2/3
    Points[1,2] = 1/3
    Points[2,1] = 1/3
    Points[2,2] = -2/3
    Points[3,1] = -2/3
    Points[3,2] = -2/3
  elseif n == 3
    NumQuad = 4 
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = -27/96
    Weights[2] = 25/96
    Weights[3] = 25/96
    Weights[4] = 25/96
    Points[1,1] = 1/3
    Points[1,2] = 1/3
    Points[2,1] = -3/5
    Points[2,2] = -3/5
    Points[3,1] = 3/5
    Points[3,2] = -3/5
    Points[4,1] = -3/5
    Points[4,2] = 3/5
  end

  return NumQuad, Weights, Points
end
