function NumberingFemDGQuad(Grid,PolyOrd)

NumGlobF = Grid.NumFaces*(PolyOrd+1)*(PolyOrd+1)

Glob = zeros(Int,(PolyOrd+1)*(PolyOrd+1),Grid.NumFaces)
GlobLoc = zeros(Int,(PolyOrd+1)*(PolyOrd+1))
for iF = 1:Grid.NumFaces
  Face = Grid.Faces[iF]
  ii = 1
  for j = 1:PolyOrd+1
    for i = 1:PolyOrd+1
      GlobLoc[ii] = i + (j - 1) * (PolyOrd + 1) + (iF - 1) * (PolyOrd + 1) *(PolyOrd + 1)
      ii = ii + 1
    end
  end
  Glob[:,iF] = GlobLoc
end

NumG = NumGlobF
NumI = NumG
Stencil = zeros(Int,Grid.NumFaces,12)

for iF = 1:Grid.NumFaces
  Stencil[iF,:] .= iF
  StencilLoc = zeros(Int, 16,1)
  StencilLoc[:] .= iF
  iS = 0
  for i = 1:4
    iN = Grid.Faces[iF].N[i]
    for j = 1:size(Grid.Nodes[iN].F,1)
      jF = Grid.Nodes[iN].F[j]
      inside = false
      for jS = 1:iS
        if StencilLoc[jS] == jF
          inside = true
          break
        end
      end
      if !inside
        iS = iS + 1
        StencilLoc[iS] = jF
      end
    end
  end
  Stencil[iF,1:iS] = StencilLoc[1:iS]
end

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



