function CurlNode!(c,uN,Metric,Grid)

  FT = eltype(c)
  t = zeros(FT,3)
  for iE = 1 : Grid.NumEdges
    N1 =  Grid.Edges[iE].N[1]
    c[N1] += Metric.DualEdge[iE]*uN[iE] 
    N2 =  Grid.Edges[iE].N[2]
    c[N2] -= Metric.DualEdge[iE]*uN[iE]
  end  
  @. c /= Metric.DualVolume
  @show maximum(Metric.DualEdge)
  @show maximum(Metric.DualVolume)
end  

function Curl1!(c,uN,Metric,Grid)

  FT = eltype(c)
  t = zeros(FT,3)
  for iF = 1 : Grid.NumFaces
    c[iF] = 0  
    for i = 1 : length(Grid.Faces[iF].E)  
      if i == 1  
        im1 = length(Grid.Faces[iF].E)
      else
        im1 = i - 1  
      end  
      if i == length(Grid.Faces[iF].E)
        ip1 = 1
      else
        ip1 = i + 1  
      end  
      iEm1 = Grid.Faces[iF].E[im1]
      uNm1 = uN[iEm1]
      nm1 = Metric.PrimalNormal[:,iEm1]
      iEp1 = Grid.Faces[iF].E[ip1]
      uNp1 = uN[iEp1]
      np1 = Metric.PrimalNormal[:,iEp1]
      iE = Grid.Faces[iF].E[i]
      t[1] = Grid.Edges[iE].t.x
      t[2] = Grid.Edges[iE].t.y
      t[3] = Grid.Edges[iE].t.z
      l = Metric.PrimalEdge[iE]
      OrientE = Grid.Faces[iF].OrientE[i]
      c[iF] += 0.5 * (uNm1 * (t' * nm1) + uNp1 * (t' * np1)) * l #* OrientE
    end    
    c[iF] /= Metric.PrimalVolume[iF]
  end
end  

function Curl!(c,uN,Metric,Grid)

  Rad = Grid.Rad
  FT = eltype(c)
  t = zeros(FT,3)
  @. c = 0
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      if i == 1
        im1 = length(Grid.Faces[iF].N)
        im2 = length(Grid.Faces[iF].N)-1
      elseif i == 2
        im1 = i - 1
        im2 = length(Grid.Faces[iF].N)
      else
        im1 = i - 1
        im2 = i - 2
      end
      if i == length(Grid.Faces[iF].N)
        ip1 = 1
      else
        ip1 = i + 1  
      end  

      iEm1 = Grid.Faces[iF].E[im1]
      uuNm1 = uN[iEm1]
      nm1 = Grid.Edges[iEm1].n
      iEm2 = Grid.Faces[iF].E[im2]
      uuNm2 = uN[iEm2]
      nm2 = Grid.Edges[iEm2].n
      iE = Grid.Faces[iF].E[i]
      uuN = uN[iE]
      n = Grid.Edges[iE].n
      iEp1 = Grid.Faces[iF].E[ip1]
      uuNp1 = uN[iEp1]
      np1 = Grid.Edges[iEp1].n
      t = Grid.Faces[iF].Mid - Grid.Edges[iE].Mid
      t = t / norm(t) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iE].Mid) * Grid.Rad
      c[iN] += (Grids.dot(t,n) * uuN + 0.5 * Grids.dot(t,nm1) * uuNm1 + 0.5 * Grids.dot(t,np1) * uuNp1)  
      t = Grid.Edges[iEm1].Mid - Grid.Faces[iF].Mid 
      t = t / norm(t) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iEm1].Mid) * Grid.Rad
      c[iN] += (0.5 * Grids.dot(t,n) * uuN + Grids.dot(t,nm1) * uuNm1 + 0.5 * Grids.dot(t,nm2) * uuNm2) 
    end
  end
  @. c /= Metric.DualVolume
end
