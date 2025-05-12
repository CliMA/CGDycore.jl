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

function NumberingFemDGTri(Grid,DoF,DoFE,PosDoFE,Proc)

  Glob = zeros(Int,DoF,Grid.NumFaces)
  GlobLoc = zeros(Int,DoF)
  for iF = 1:Grid.NumFaces
    Face = Grid.Faces[iF]
    for iDoF = 1 : DoF
      for i = 1 : DoF
        GlobLoc[iDoF] = iDoF + (iF - 1) * DoF
      end
    end
    @. Glob[:,iF] = GlobLoc
  end

  GlobE = zeros(Int,2,DoFE,Grid.NumEdges)
  GlobLocE = zeros(Int,2,DoFE)
  iEB = 0
  for iE = 1 : Grid.NumEdges
    iF = Grid.Edges[iE].F[1]  
    if iF > 0 
      FE = Grid.Edges[iE].FE[1]
      if iF <= Grid.NumFaces 
        for iDoFE = 1 : DoFE
          GlobLocE[1,iDoFE] = PosDoFE[iDoFE,FE] + (iF - 1) * DoF
        end
      else
        iEB += 1
        for iDoFE = 1 : DoFE
          GlobLocE[1,iDoFE] = iDoFE + (iEB - 1) * DoFE + Grid.NumFaces * DoF
        end
      end  
    end  
    iF = Grid.Edges[iE].F[2]  
    if iF > 0 
      FE = Grid.Edges[iE].FE[2]
      if iF <= Grid.NumFaces 
        for iDoFE = 1 : DoFE
          GlobLocE[2,iDoFE] = PosDoFE[DoFE-iDoFE+1,FE] + (iF - 1) * DoF
        end
      else  
        iEB += 1
        for iDoFE = 1 : DoFE
          GlobLocE[2,iDoFE] = DoFE - iDoFE + 1 + (iEB - 1) * DoFE + Grid.NumFaces * DoF
        end
      end  
    end  
    @. GlobE[:,:,iE] = GlobLocE
  end  

  NumI = Grid.NumFaces * DoF
  NumG = NumI + Grid.NumEdgesB * DoFE
  Stencil = zeros(Int,0,0)

  MasterSlave = zeros(Int,NumG) # Sinn
  BoundaryDoF = zeros(Int,0)
    
  return (Glob,GlobE,NumG,NumI,Stencil,MasterSlave,BoundaryDoF)
end



