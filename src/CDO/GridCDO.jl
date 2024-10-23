function InitGridCDO(FT,Grid)
  nz = Grid.nz
  NumNodes = Grid.NumNodes * (nz + 1)
  NumEdges = Grid.NumEdges * (nz + 1) + Grid.NumNodes * nz
  NumFaces = Grid.NumFaces * (nz + 1) + Grid.NumEdges * nz
  NumCells = Grid.NumFaces * nz
  # Node numbering
  # Nodes in the first layer, second layer and so on
  # Edge numbering 
  # Edges in the first layer, second layer and so on
  # Nodes in the first layer, second layer and so on
  # Face numbering 
  # Faces in the first layer, second layer and so on
  # Edges in the first layer, second layer and so on
  # Cell numbering 
  # Faces in the first layer, second layer and so on

  # Incidence
  # Edge Node
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Int32[]
  for iEG = 1 : Grid.NumEdges 
    iNG1 = Grid.Edges[iEG].N[1]  
    iNG2 = Grid.Edges[iEG].N[2]  
    for iz = 1 : nz + 1
      iE = iEG + (iz - 1) * Grid.NumEdges  
      iN1 = iNG1 + (iz - 1) * Grid.NumNodes
      iN2 = iNG2 + (iz - 1) * Grid.NumNodes
      push!(RowInd,iE)
      push!(ColInd,iN1)
      push!(Val,-1)
      push!(RowInd,iE)
      push!(ColInd,iN2)
      push!(Val,1)
    end  
  end  
  for iNG = 1 : Grid.NumNodes 
    for iz = 1 : nz
      iE = iNG + (iz - 1) *Grid.NumNodes + Grid.NumEdges * (nz + 1)   
      iN1 = iNG + (iz - 1) * Grid.NumNodes 
      iN2 = iNG + iz * Grid.NumNodes 
      push!(RowInd,iE)
      push!(ColInd,iN1)
      push!(Val,-1)
      push!(RowInd,iE)
      push!(ColInd,iN2)
      push!(Val,1)
    end  
  end  
  IncEN = sparse(RowInd, ColInd, Val)
  # Incidence
  # Face Edge  
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Int32[]
  for iFG = 1 : Grid.NumFaces 
    for iz = 1 : nz + 1
      iF = iFG + (iz - 1) * Grid.NumFaces  
      for i = 1 : length(Grid.Faces[iFG].E)
        iEG = Grid.Faces[iFG].E[i]
        iE = iEG + (iz - 1) * Grid.NumEdges  
        push!(RowInd,iF)
        push!(ColInd,iE)
        push!(Val,Grid.Faces[iFG].OrientE[i])
      end
    end
  end  
  for iEG = 1 : Grid.NumEdges
    iNG1 = Grid.Edges[iEG].N[1]  
    iNG2 = Grid.Edges[iEG].N[2]  
    for iz = 1 : nz
      iF = iEG + (iz - 1) * Grid.NumEdges + Grid.NumFaces * (nz + 1) 
      iE = iEG + (iz - 1) * Grid.NumEdges 
      push!(RowInd,iF) 
      push!(ColInd,iE)
      push!(Val,1)
      iE = iNG2 + (iz - 1) * Grid.NumNodes + (nz + 1) * Grid.NumEdges
      push!(RowInd,iF) 
      push!(ColInd,iE)
      push!(Val,1)
      iE = iEG + iz * Grid.NumEdges 
      push!(RowInd,iF) 
      push!(ColInd,iE)
      push!(Val,-1)
      iE = iNG1 + (iz - 1) * Grid.NumNodes + (nz + 1) * Grid.NumEdges
      push!(RowInd,iF) 
      push!(ColInd,iE)
      push!(Val,-1)
    end
  end  
  IncFE = sparse(RowInd, ColInd, Val)
  # Incidence
  # Cell Face  
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Int32[]
  for iFG = 1 : Grid.NumFaces
    for iz = 1 : nz
      iC = iFG + (iz - 1) * Grid.NumFaces
      iF = iFG + (iz - 1) * Grid.NumFaces
      push!(RowInd,iC) 
      push!(ColInd,iF)
      push!(Val,-1)
      iF = iFG + iz * Grid.NumFaces
      push!(RowInd,iC) 
      push!(ColInd,iF)
      push!(Val,1)
      for i = 1 : length(Grid.Faces[iFG].E)
        iEG = Grid.Faces[iFG].E[i]  
        iF = iEG + (iz - 1) * Grid.NumEdges + (nz + 1) * Grid.NumFaces
        push!(RowInd,iC) 
        push!(ColInd,iF)
        push!(Val,Grid.Faces[iFG].OrientE[i])
      end
    end  
  end  
  IncCF = sparse(RowInd, ColInd, Val)

  return GridCDO{FT,
                 typeof(IncEN)}(
    IncEN,
    IncFE,
    IncCF,
    )
end  

        



