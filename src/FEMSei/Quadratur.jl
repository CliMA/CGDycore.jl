function QuadRule(type::Grids.QuadPrimal,n)
  NumQuad = n * n
  x, w = gaussradau(n)
  Weights = zeros(NumQuad)
  Points = zeros(NumQuad,2)
  ii = 1
  for j = 1 : n
    for i = 1 : n
      Weights[ii] = w[i] * w[j]
      Points[ii,1] = x[i]
      Points[ii,2] = x[j]
      ii += 1
    end
  end
  return NumQuad, Weights, Points
end


function QuadRule(type::Grids.QuadDual,n)
  NumQuad =  n * n  
  x, w = gaussradau(n)
  x .= -x
  Weights = zeros(NumQuad)
  Points = zeros(NumQuad,2)
  ii = 1
  for j = 1 : n
    for i = 1 : n
      Weights[ii] = w[n+1-i] * w[n+1-j]
      Points[ii,1] = x[i]
      Points[ii,2] = x[j]
      ii += 1
    end
  end
  return NumQuad, Weights, Points
end

function QuadRule(type::Grids.Line,n)
 x, w = gausslobatto(n+1)
# x, w = gausslegendre(n+1)
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
# x, w = gausslegendre(n+1)
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

function QuadRule(type::Grids.Tri,Ord)

  if Ord == 1
    NumQuad = 1  
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = 1
    Points[1,1] = 1/3
    Points[1,2] = 1/3
    @. Weights = 2.0 * Weights
    for i = 1 : NumQuad
      Points[i,1] = 2.0 * Points[i,1] - 1.0  
      Points[i,2] = 2.0 * Points[i,2] - 1.0  
    end  

  elseif Ord == 2
    NumQuad = 3  
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = 1/3
    Weights[2] = 1/3
    Weights[3] = 1/3
    Points[1,1] = 1/6
    Points[1,2] = 2/3
    Points[2,1] = 1/6
    Points[2,2] = 1/6
    Points[3,1] = 2/3
    Points[3,2] = 1/6
    @. Weights = 2.0 * Weights
    for i = 1 : NumQuad
      Points[i,1] = 2.0 * Points[i,1] - 1.0  
      Points[i,2] = 2.0 * Points[i,2] - 1.0  
    end  
  
  elseif Ord == 3
    NumQuad = 4 
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = -27/96
    Weights[2] = 25/96
    Weights[3] = 25/96
    Weights[4] = 25/96
    Points[1,1] = 1/3
    Points[1,2] = 1/3
    Points[2,1] = 0.2
    Points[2,2] = 0.2
    Points[3,1] = 0.6
    Points[3,2] = 0.2
    Points[4,1] = 0.2
    Points[4,2] = 0.6
    @. Weights = 4.0 * Weights
    for i = 1 : NumQuad
      Points[i,1] = 2.0 * Points[i,1] - 1.0  
      Points[i,2] = 2.0 * Points[i,2] - 1.0  
    end  
  elseif Ord == 4
    NumQuad = 6
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = 0.16666666666666666667
    Weights[2] = 0.16666666666666666667
    Weights[3] = 0.16666666666666666667
    Weights[4] = 0.16666666666666666667
    Weights[5] = 0.16666666666666666667
    Weights[6] = 0.16666666666666666667
    Points[1,1] = 0.659027622374092  
    Points[1,2] = 0.231933368553031
    Points[2,1] = 0.659027622374092  
    Points[2,2] = 0.109039009072877
    Points[3,1] = 0.231933368553031  
    Points[3,2] = 0.659027622374092
    Points[4,1] = 0.231933368553031  
    Points[4,2] = 0.109039009072877
    Points[5,1] = 0.109039009072877  
    Points[5,2] = 0.659027622374092
    Points[6,1] = 0.109039009072877  
    Points[6,2] = 0.231933368553031
    @. Weights = 4.0 * Weights
    for i = 1 : NumQuad
      Points[i,1] = 2.0 * Points[i,1] - 1.0  
      Points[i,2] = 2.0 * Points[i,2] - 1.0  
    end  
  elseif Ord == 5
    NumQuad = 6
    Weights = zeros(NumQuad)
    Points = zeros(NumQuad,2)
    Weights[1] = 0.109951743655322
    Weights[2] = 0.109951743655322
    Weights[3] = 0.109951743655322
    Weights[4] = 0.223381589678011
    Weights[5] = 0.223381589678011
    Weights[6] = 0.223381589678011
    Points[1,1] = 0.816847572980459 
    Points[1,2] = 0.091576213509771
    Points[2,1] = 0.091576213509771 
    Points[2,2] = 0.816847572980459
    Points[3,1] = 0.091576213509771 
    Points[3,2] = 0.091576213509771
    Points[4,1] = 0.108103018168070 
    Points[4,2] = 0.445948490915965
    Points[5,1] = 0.445948490915965 
    Points[5,2] = 0.108103018168070
    Points[6,1] = 0.445948490915965 
    Points[6,2] = 0.445948490915965
    @. Weights = 2.0 * Weights
    for i = 1 : NumQuad
      Points[i,1] = 2.0 * Points[i,1] - 1.0  
      Points[i,2] = 2.0 * Points[i,2] - 1.0  
    end  
  end

  return NumQuad, Weights, Points
end
