function Gradient(backend,FTB,MetricFV,Grid)
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]

  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      push!(ColInd,Grid.Faces[iF].F)
      iE =  Grid.Faces[iF].E[i]
      push!(RowInd,iE)
      Value = Grid.Faces[iF].OrientE[i]
      push!(Val,Value)
    end  
  end
  Grad = sparse(RowInd, ColInd, Val)
  return Grad
end
