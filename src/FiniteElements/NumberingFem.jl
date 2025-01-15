function NumberingFemCGTri(Grid,PolyOrd,DoF)
  if PolyOrd==1
#   1
#   |      
#   2---3   
   Per=[1 2 3];
  end    
  Glob=zeros(Int,Dof,Grid.NumFaces)
  GlobLoc=zeros(Int,DoF)
  for iF=1:Grid.NumFaces
    Face=Grid.Faces[iF];
    ii=1;
    for i=1:size(Face.N,1)
      GlobLoc[ii]=Grid.Nodes[Face.N[i]].N;
      ii=ii+1;
    end
    Glob[:,iF]=GlobLoc
  end
end

function NumberingFemCGQuad(Grid,PolyOrd)
  if PolyOrd==0
#   1
    Per=[1];
  elseif PolyOrd==1
#    4    3
#    1    2
    Per=[1 2 4 3];
  elseif PolyOrd==2
#   4   7   3
#   8   9   6
#   1   5   2
    Per=[1 5 2 8 9 6 4 7 3];
  elseif PolyOrd==3
#    4 10   9 3
#   11 15  16 8
#   12 13  14 7
#    1  5   6  2
    Per=[1 5 6 2 12 13 14 7 11 15 16 8 4 10 9 3];
  elseif PolyOrd==4
#    4 13 12 11  3
#   14 23 24 25 10
#   15 20 21 22  9
#   16 17 18 19  8
#    1  5  6  7  2
    Per=[1  5  6  7  2 16 17 18 19  8 15 20 21 22  9 14 23 24 25 10 4 13 12 11  3];
  else
    nIn=PolyOrd-1;
    Row1=5:1:5+nIn-1;
    Col2=5+2*nIn-1:-1:5+nIn;
    Row2=5+3*nIn-1:-1:5+2*nIn;
    Col1=5+3*nIn:1:5+4*nIn-1;
    Core=zeros(nIn,nIn);
    Ind=5+4*nIn;
    for i=nIn:-1:1
      for j=1:nIn
        Core(i,j)=Ind;
        Ind=Ind+1;
      end
    end
    Mat=[  4   Row2  3
         Col1' Core Col2'
           1   Row1   2];
    Per=[];
    for i=PolyOrd+1:-1:1
      Per=[Per Mat(i,:)];
    end
  end
  Per = vec(Per);
  NumGlobN=Grid.NumNodes;
  NumGlobE=Grid.NumEdges*(PolyOrd-1);
  NumGlobF=Grid.NumFaces*(PolyOrd-1)*(PolyOrd-1);

  Glob=zeros(Int,(PolyOrd+1)*(PolyOrd+1),Grid.NumFaces)
  GlobLoc=zeros(Int,(PolyOrd+1)*(PolyOrd+1))
  for iF=1:Grid.NumFaces
    Face=Grid.Faces[iF];
    ii=1;
    for i=1:size(Face.N,1)
      GlobLoc[ii]=Grid.Nodes[Face.N[i]].N;
      ii=ii+1;
    end
    for i=1:size(Face.E,1)
      iE=Face.E[i];
      if Grid.Edges[iE].N[1]==Face.N[i]
        for j=1:PolyOrd-1
          GlobLoc[ii]=(Grid.Edges[Face.E[i]].E-1)*(PolyOrd-1)+j+NumGlobN;
          ii=ii+1;
        end
      else
        for j=PolyOrd-1:-1:1
          GlobLoc[ii]=(Grid.Edges[Face.E[i]].E-1)*(PolyOrd-1)+j+NumGlobN;
          ii=ii+1;
        end
      end
    end
    for j=1:PolyOrd-1
      for i=1:PolyOrd-1
        GlobLoc[ii]=NumGlobN+NumGlobE+i+(j-1)*(PolyOrd-1)+(iF-1)*(PolyOrd-1)*(PolyOrd-1);
        ii=ii+1;
      end
    end
    GlobLoc=GlobLoc[Per]
    Glob[:,iF]=GlobLoc
  end
  # Boundary list 
  BoundaryDoF = zeros(Int,0)
  for iN = 1 : Grid.NumNodes
    if Grid.Nodes[iN].Type == 'B'  
      push!(BoundaryDoF,Grid.Nodes[iN].N)
    end  
  end  
  for iF = 1 : Grid.NumFaces
    Face = Grid.Faces[iF]
    for i = 1 : size(Face.E,1)
      iE = Face.E[i]
      if Grid.Edges[iE].Type == "B"
        if Grid.Edges[iE].N[1]==Face.N[i] 
          for j = 1 : PolyOrd - 1
            push!(BoundaryDoF,(Grid.Edges[Face.E[i]].E-1)*(PolyOrd-1)+j+NumGlobN)
          end
        else
          for j = PolyOrd - 1 : -1 : 1
            push!(BoundaryDoF,(Grid.Edges[Face.E[i]].E-1)*(PolyOrd-1)+j+NumGlobN)
          end
        end
      end
    end
  end  

  NumG=NumGlobN+NumGlobE+NumGlobF;
  NumI=NumG;
  Stencil=zeros(Int,Grid.NumFaces,13);

  for iF=1:Grid.NumFaces
    Stencil[iF,:] .= iF;
    StencilLoc=zeros(Int, 16,1);
    StencilLoc[:] .= iF;
    iS=0;
    for i = 1 : length(Grid.Faces[iF].N)
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
    Stencil[iF,1:iS]=StencilLoc[1:iS];
  end

  MasterSlave=zeros(Int,NumG)
  ii=1
  for iN=1:Grid.NumNodes
    MasterSlave[ii] = Grid.Nodes[iN].MasterSlave
    ii = ii + 1
  end  
  for iE = 1 : Grid.NumEdges
    for i=1:PolyOrd-1
      MasterSlave[ii] = Grid.Edges[iE].MasterSlave 
      ii = ii +1
    end
  end  
  for iF = 1 : Grid.NumFaces
    for j=1:PolyOrd-1
      for i=1:PolyOrd-1
        MasterSlave[ii] = 1
        ii = ii +1
      end  
    end
  end  
    
  return (Glob,NumG,NumI,Stencil,MasterSlave,BoundaryDoF)
end

function NumberingFemRT0(Grid,PolyOrd)
NumGlobN=0
NumGlobE=Grid.NumEdges
NumGlobF=0

Glob=zeros(Int,3,1,Grid.NumFaces)
  GlobLoc=zeros(Int,3)
  for iF=1:Grid.NumFaces
    Face=Grid.Faces[iF];
    ii=1;
    for i=1:size(Face.E,1)
      iE=Face.E[i];
      if Grid.Edges[iE].N[1]==Face.N[i]
        for j=1:PolyOrd-1
          GlobLoc[ii]=Grid.Edges[Face.E[i]].E
          ii=ii+1;
        end
      end
    end
  end
  Glob[:,:,iF]=reshape(GlobLoc,3,1);

  NumG=NumGlobN+NumGlobE+NumGlobF;
  NumI=NumG;
  Stencil=zeros(Int,Grid.NumFaces,13);

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
    Stencil[iF,1:iS]=StencilLoc[1:iS];
  end

  MasterSlave=zeros(Int,NumG)
  ii=1
  for iN=1:Grid.NumNodes
    MasterSlave[ii] = Grid.Nodes[iN].MasterSlave
    ii = ii + 1
  end  
  for iE = 1 : Grid.NumEdges
    for i=1:PolyOrd-1
      MasterSlave[ii] = Grid.Edges[iE].MasterSlave 
      ii = ii +1
    end
  end  
  for iF = 1 : Grid.NumFaces
    for j=1:PolyOrd-1
      for i=1:PolyOrd-1
        MasterSlave[ii] = 1
        ii = ii +1
      end  
    end
  end  
    
  return (Glob,NumG,NumI,Stencil,MasterSlave)
end

