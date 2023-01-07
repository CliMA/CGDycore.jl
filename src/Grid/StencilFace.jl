function StencilFace(Grid)
  Stencil=zeros(Int,Grid.NumFaces,12);

  for iF=1:Grid.NumFaces
    Stencil[iF,:] .= iF;
    StencilLoc=zeros(Int, 16,1);
    StencilLoc[:] .= iF;
    iS=0;
    for i=1:4
      iN=Grid.Faces[iF].N[i];
      for j=1:size(Grid.Nodes[iN].F,1)
        jF=Grid.Nodes[iN].F[j];
        inside=false;
        for jS=1:iS
          if StencilLoc[jS]==jF
            inside=true;
            break
          end
        end
        if !inside
          iS=iS+1;
          StencilLoc[iS]=jF;
        end
      end
    end
    Grid.Faces[iF].Stencil = zeros(Int,iS)
    Grid.Faces[iF].Stencil .= StencilLoc[1:iS]
  end
  return Grid
end
