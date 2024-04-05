function Divergence(backend,FTB,MetricFV,Grid)
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]

  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      push!(RowInd,Grid.Faces[iF].F)
      iE =  Grid.Faces[iF].E[i]
      push!(ColInd,iE)
      Value = Grid.Faces[iF].OrientE[i] * MetricFV.PrimalEdge[iE] / MetricFV.PrimalVolume[iF]
      push!(Val,Value)
    end  
  end
  Div = sparse(RowInd, ColInd, Val)
  return Div
end
