function NumberingFemDGQuad(Grid,PolyOrd,Proc)

  N = PolyOrd + 1
  Glob = zeros(Int,N*N,Grid.NumFaces)
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
    @. Glob[:,iF] = GlobLoc
  end

  GlobE = zeros(Int,2,N,Grid.NumEdges)
  GlobLocE = zeros(Int,2,N)
  iEB = 0
  for iE = 1 : Grid.NumEdges
    for j = 1 : 2  
      iF = Grid.Edges[iE].F[j]  
      if iF > 0 
        FE = Grid.Edges[iE].FE[j]
        if iF <= Grid.NumFaces 
          if FE == 1  
            ii = 1  
            for i = 1 : N
              GlobLocE[j,i] = ii + (iF - 1) * N * N
              ii += 1
            end
          elseif FE == 2  
            ii =  N
            for i = 1 : N
              GlobLocE[j,i] = ii + (iF - 1) * N * N
              ii += N
            end
          elseif FE == 3  
            ii =  1 + (N - 1) * N 
            for i = 1 : N
              GlobLocE[j,i] = ii + (iF - 1) * N * N
              ii += 1
            end  
          elseif FE == 4  
            ii =  1 
            for i = 1 : N
              GlobLocE[j,i] = ii + (iF - 1) * N * N
              ii += N
            end  
          end
        else
          iEB += 1
          for i = 1 : N
            GlobLocE[j,i] = i + (iEB - 1) * N + Grid.NumFaces * N * N
          end
        end
      end  
    end  
    @. GlobE[:,:,iE] = GlobLocE
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
    
  return (Glob,GlobE,NumG,NumI,Stencil,MasterSlave,BoundaryDoF)
end



