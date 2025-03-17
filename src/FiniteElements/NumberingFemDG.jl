function NumberingFemDGQuad(Grid,PolyOrd)

  N = PolyOrd + 1
  Glob = zeros(Int,N*N,Grid.NumFaces+Grid.NumEdgesB)
  GlobLoc = zeros(Int,N*N)
  for iF = 1:Grid.NumFaces
    Face = Grid.Faces[iF]
    ii = 1
    for j = 1 : N
      for i = 1 : N
        GlobLoc[ii] = i + (j - 1) * N + (iF - 1) * N * N
        ii = ii + 1
      end
    end
    Glob[:,iF] = GlobLoc
  end
  for iE = 1 : Grid.NumEdgesB
    @. GlobLoc = 0  
    if Grid.Edges[iE].F[1] < Grid.Edges[iE].F[2]
      iF = Grid.Edges[iE].F[2]
      FE = Grid.Edges[iE].FE[2]
    else  
      iF = Grid.Edges[iE].F[1]
      FE = Grid.Edges[iE].FE[1]
    end  
    if FE == 1
      ii = 1  
      for i = 1 : N
        GlobLoc[ii] = i + (iE - 1) * N + Grid.NumFaces * N * N
        ii += 1
      end
    elseif FE == 2  
      ii =  N
      for i = 1 : PolyOrd + 1
        GlobLoc[ii] = i + (iE - 1) * N + Grid.NumFaces * N * N
        ii += N
      end
    elseif FE == 3  
      ii =  1 + PolyOrd * N 
      for i = 1 : N
        GlobLoc[ii] = i + (iE - 1) * N + Grid.NumFaces * N * N
        ii += 1
      end  
    elseif FE == 4  
      ii =  1 
      for i = 1 : PolyOrd + 1
        GlobLoc[ii] = i + (iE - 1) * N + Grid.NumFaces * N * N
        ii += PolyOrd + 1
      end  
    end
    Glob[:,iF] = GlobLoc
  end  

  NumI = Grid.NumFaces * N * N
  NumG = NumI + Grid.NumEdgesB * N
  Stencil = zeros(Int,0,0)

  MasterSlave = zeros(Int,NumG)
  ii = 1
  for iF = 1 : Grid.NumFaces
    for j = 1 : PolyOrd - 1
      for i = 1 : PolyOrd - 1
        MasterSlave[ii] = 1
        ii = ii + 1
      end  
    end
  end  
  BoundaryDoF = zeros(Int,0)
    
  return (Glob,NumG,NumI,Stencil,MasterSlave,BoundaryDoF)
end



