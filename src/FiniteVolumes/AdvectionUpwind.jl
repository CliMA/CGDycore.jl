function AdvectionUpwind!(backend,FTB,Rhs,p,u,MetricFV,Grid)

  for iF = 1 : Grid.NumFaces
#   @show Grid.Faces[iF].Orientation  
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
#     @show iE,Grid.Edges[iE].F
      if u[iE]*Grid.Faces[iF].OrientE[i] > 0.0
        Rhs[iF] -= u[iE] * Grid.Faces[iF].OrientE[i] * p[iF] *
          MetricFV.PrimalEdge[iE] / MetricFV.PrimalVolume[iF]   
        if iF == Grid.Edges[iE].F[1]
          iFN = Grid.Edges[iE].F[2]
        else 
          iFN = Grid.Edges[iE].F[1]  
        end  
        Rhs[iFN] += u[iE] * Grid.Faces[iF].OrientE[i] * p[iF] *
          MetricFV.PrimalEdge[iE] / MetricFV.PrimalVolume[iFN]
      end
    end
  end
end

function AdvectionUpwindE!(backend,FTB,Rhs,p,u,MetricFV,Grid)

  list = [10000 10011 9997 9986]

  for iE = 1 : Grid.NumEdges
    if length(Grid.Edges[iE].F) > 1  
      iFL = Grid.Edges[iE].F[1]  
      iFR = Grid.Edges[iE].F[2]  
      MomFlux = u[iE] * MetricFV.PrimalEdge[iE]
      if MomFlux > 0  
        Rhs[iFL] -= MomFlux * p[iFL]  
        Rhs[iFR] += MomFlux * p[iFL]  
      else
        Rhs[iFL] -= MomFlux * p[iFR]  
        Rhs[iFR] += MomFlux * p[iFR]  
      end
    end
  end
  @. Rhs /= MetricFV.PrimalVolume
end

