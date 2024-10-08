function DivMPFA(backend,FT,Metric,Grid)
  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]

  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i] 
      push!(RowInd,Grid.Faces[iF].F)
      push!(ColInd,Grid.Edges[iE].E)
      push!(Val,Grid.Faces[iF].OrientE[i] * Metric.PrimalEdge[iE] / Metric.PrimalVolume[iF])
    end  
  end    
  Div = sparse(RowInd, ColInd, Val)
end

function GradMPFA1(KiteGrid,Grid)

  Nodes = KiteGrid.Nodes
  Faces = KiteGrid.Faces
  NumNodes = KiteGrid.NumNodes

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  Shift = Grid.NumNodes + Grid.NumEdges
  for iN = 1 : NumNodes
    if Nodes[iN].Type == 'N'  
      NumF = length(Nodes[iN].F) 
      I = zeros(2 * NumF +1,2 * NumF +1)
      ILoc = zeros(4,4,NumF)
      GlobLoc = zeros(Int,4,NumF)
      b = zeros(2 * NumF + 1)
      c = zeros(2 * NumF + 1)
      Grad = zeros(NumF,NumF)
      for i = 1 : NumF
        GlobLoc[:,i] = [i,NumF + i,2 * NumF + 1,NumF + i - 1]  
        if i == 1
          GlobLoc[4,i] = 2 * NumF  
        end  
        iF = Nodes[iN].F[i]  
        @views LocalInterpolation1!(ILoc[:,:,i],Faces[iF],KiteGrid)
        @views @. I[GlobLoc[:,i],GlobLoc[:,i]] += ILoc[:,:,i]
      end
      @. I[end,1:NumF] = 1
      I[end,end] = -NumF
      for i = 1 : NumF
        iF = Nodes[iN].F[i]  
        @. b = 0
        b[i] = 1
        c = I \ b
        for j = 1 : NumF
          jF = Nodes[iN].F[j]  
          Temp = sum(ILoc[2,:,j] .* c[GlobLoc[:,j]])
          push!(RowInd,div(KiteGrid.Faces[jF].E[2]+1,2))
          push!(ColInd,KiteGrid.Faces[iF].N[1]-Shift)
          push!(Val,Temp)
        end  
      end  
    end  
  end  
  Grad = sparse(RowInd, ColInd, Val)
end

function GradMPFA(backend,FT,Grid)

  Nodes = Grid.Nodes
  Edges = Grid.Edges
  Faces = Grid.Faces
  NumNodes = Grid.NumNodes
  NumFaces = Grid.NumFaces

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  ValGrad = Float64[]
  Shift = Grid.NumNodes + Grid.NumEdges
  for iN = 1 : NumNodes
    NG = Grid.Nodes[iN].NG  
    EdgesInNode = Int64[]
    NumF = length(Nodes[iN].F) 
    I = zeros(2 * NumF +1,2 * NumF +1)
    ILoc = zeros(4,4,NumF)
    GlobLoc = zeros(Int,4,NumF)
    b = zeros(2 * NumF + 1)
    c = zeros(2 * NumF + 1)
    for i = 1 : NumF
      iF = Nodes[iN].F[i]
      iE2 = 0
      iE3 = 0
      # Find Node in Face
      for j = 1 : length(Faces[iF].N)
        if iN == Faces[iF].N[j]
          if j == 1
            iE2 = Faces[iF].E[end]
            push!(EdgesInNode,length(Faces[iF].N))
          else
            iE2 = Faces[iF].E[j-1]
            push!(EdgesInNode,j-1)
          end  
          iE3 =  Faces[iF].E[j]
          exit
        end  
      end
      GlobLoc[:,i] = [i,NumF + i,2 * NumF + 1,NumF + i - 1]  
      if i == 1
        GlobLoc[4,i] = 2 * NumF  
      end  
      P1 = Faces[iF].Mid
      P2 = Edges[iE2].Mid
      P3 = Nodes[iN].P
      P4 = Edges[iE3].Mid
      @views LocalInterpolationGrad!(ILoc[:,:,i],P1,P2,P3,P4,Grid.Rad)
      @views @. I[GlobLoc[:,i],GlobLoc[:,i]] += ILoc[:,:,i]
    end
    @. I[end,1:NumF] = 1
    I[end,end] = -NumF
    for i = 1 : NumF
      iF = Nodes[iN].F[i]  
      @. b = 0
      b[i] = 1
      c = I \ b
      for j = 1 : NumF
        jF = Nodes[iN].F[j]  
        FG = Faces[iF].FG
        jE = Faces[jF].E[EdgesInNode[j]]
        EG = Edges[jE].EG
        TempGrad = -sum(ILoc[2,:,j] .* c[GlobLoc[:,j]]) * Faces[jF].OrientE[EdgesInNode[j]]
#       TempGrad = -sum(ILoc[2,:,j] .* c[GlobLoc[:,j]]) 
        Temp = 0.5 * c[GlobLoc[2,j]] 
        push!(RowInd,jE)
        push!(ColInd,iF)
        push!(ValGrad,TempGrad)
        push!(Val,Temp)
      end  
    end  
  end  
  Grad = sparse(RowInd, ColInd, ValGrad)
  Inter = sparse(RowInd, ColInd, Val)
  return Grad[1:Grid.NumEdges,:],Inter[1:Grid.NumEdges,:]
end
function GradMPFATri(backend,FT,Grid)

  Nodes = Grid.Nodes
  Edges = Grid.Edges
  Faces = Grid.Faces
  NumNodes = Grid.NumNodes
  NumFaces = Grid.NumFaces

  RowInd = Int64[]
  ColInd = Int64[]
  Val = Float64[]
  ValGrad = Float64[]
  Shift = Grid.NumNodes + Grid.NumEdges
  n1 = zeros(3)
  n2 = zeros(3)
  for iN = 1 : NumNodes
    NG = Grid.Nodes[iN].NG  
    EdgesInNode = Int64[]
    NumF = length(Nodes[iN].F) 
    I = zeros(2 * NumF,2 * NumF)
    ILoc = zeros(3,3,NumF)
    GlobLoc = zeros(Int,3,NumF)
    b = zeros(2 * NumF)
    c = zeros(2 * NumF)
    for i = 1 : NumF
      iF = Nodes[iN].F[i]
      iE2 = 0
      iE3 = 0
      # Find Node in Face
      for j = 1 : length(Faces[iF].N)
        if iN == Faces[iF].N[j]
          if j == 1
            iE2 = Faces[iF].E[end]
            push!(EdgesInNode,length(Faces[iF].N))
          else
            iE2 = Faces[iF].E[j-1]
            push!(EdgesInNode,j-1)
          end  
          iE3 =  Faces[iF].E[j]
          exit
        end  
      end
      GlobLoc[:,i] = [i,NumF + i,NumF + i - 1]  
      if i == 1
        GlobLoc[3,i] = 2 * NumF  
      end  
      P1 = Faces[iF].Mid
      P2 = Edges[iE2].Mid
      n1[1] = Edges[iE2].n.x
      n1[2] = Edges[iE2].n.y
      n1[3] = Edges[iE2].n.z
#     P3 = Nodes[iN].P
      P3 = Edges[iE3].Mid
      n2[1] = Edges[iE3].n.x
      n2[2] = Edges[iE3].n.y
      n2[3] = Edges[iE3].n.z
      @views LocalInterpolationGradTri!(ILoc[:,:,i],P1,P2,P3,n1,n2,Grid.Rad)
      @views @. I[GlobLoc[:,i],GlobLoc[:,i]] += ILoc[:,:,i]
    end
    for i = 1 : NumF
      iF = Nodes[iN].F[i]  
      @. b = 0
      b[i] = 1
      c = I \ b
      for j = 1 : NumF
        jF = Nodes[iN].F[j]  
        FG = Faces[iF].FG
        jE = Faces[jF].E[EdgesInNode[j]]
        EG = Edges[jE].EG
        TempGrad = 0.5 * sum(ILoc[2,:,j] .* c[GlobLoc[:,j]])  #* Faces[jF].OrientE[EdgesInNode[j]]
        Temp = 0.5 * c[GlobLoc[2,j]] 
        push!(RowInd,jE)
        push!(ColInd,iF)
        push!(ValGrad,TempGrad)
        push!(Val,Temp)
      end  
    end  
  end  
  Grad = sparse(RowInd, ColInd, ValGrad)
  Inter = sparse(RowInd, ColInd, Val)
  return Grad[1:Grid.NumEdges,:],Inter[1:Grid.NumEdges,:]
end

function LocalInterpolation1!(ILoc,Face,Grid)

  grad = zeros(2)
  n = zeros(2)
  @polyvar x y
  lx0 = 0.5*(1-x)
  lx1 = 0.5*(1+x)
  ly0 = 0.5*(1-y)
  ly1 = 0.5*(1+y)
  phi1 = lx0 * ly0
  phi2 = lx1 * ly0
  phi3 = lx1 * ly1
  phi4 = lx0 * ly1
  gradphi1x = differentiate(phi1,x)
  gradphi1y = differentiate(phi1,y)
  gradphi2x = differentiate(phi2,x)
  gradphi2y = differentiate(phi2,y)
  gradphi3x = differentiate(phi3,x)
  gradphi3y = differentiate(phi3,y)
  gradphi4x = differentiate(phi4,x)
  gradphi4y = differentiate(phi4,y)
  @. ILoc = 0
  # Pointwise interpolation
  ILoc[1,1] = 1
  #First gradient

  ksi1 = 1.0
  ksi2 = -1.0
  n[1] = 1.0
  n[2] = 0.0
  P1 = Grid.Nodes[Face.N[2]].P
  P2 = Grid.Nodes[Face.N[3]].P
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  FEMSei.Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Face,Grid)
  eInv = 2.0 * J / (Grids.SizeGreatCircle(P1,P2) * Grid.Rad)
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[2,1] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[2,2] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[2,3] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi4x(ksi1,ksi2) 
  grad[2] = gradphi4y(ksi1,ksi2) 
  ILoc[2,4] = eInv * (pinvDF * n)' * (pinvDF * grad)
  ksi1 = -1.0
  ksi2 = 1.0
  n[1] = 0.0
  n[2] = 1.0
  P1 = Grid.Nodes[Face.N[3]].P
  P2 = Grid.Nodes[Face.N[4]].P
  FEMSei.Jacobi(DF,detDF,pinvDF,X,Grid.Type,ksi1,ksi2,Face,Grid)
  eInv = 2.0 * J / (Grids.SizeGreatCircle(P1,P2) * Grid.Rad)
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[4,1] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[4,2] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[4,3] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi4x(ksi1,ksi2) 
  grad[2] = gradphi4y(ksi1,ksi2) 
  ILoc[4,4] = eInv * (pinvDF * n)' * (pinvDF * grad)
end

function LocalInterpolationGrad!(ILoc,P1,P2,P3,P4,Rad)

  grad = zeros(2)
  n = zeros(2)
  @polyvar x y
  lx0 = 0.5*(1-x)
  lx1 = 0.5*(1+x)
  ly0 = 0.5*(1-y)
  ly1 = 0.5*(1+y)
  phi1 = lx0 * ly0
  phi2 = lx1 * ly0
  phi3 = lx1 * ly1
  phi4 = lx0 * ly1
  gradphi1x = differentiate(phi1,x)
  gradphi1y = differentiate(phi1,y)
  gradphi2x = differentiate(phi2,x)
  gradphi2y = differentiate(phi2,y)
  gradphi3x = differentiate(phi3,x)
  gradphi3y = differentiate(phi3,y)
  gradphi4x = differentiate(phi4,x)
  gradphi4y = differentiate(phi4,y)
  @. ILoc = 0
  # Pointwise interpolation
  ILoc[1,1] = 1
  #First gradient

  ksi1 = 1.0
  ksi2 = -1.0
  n[1] = 1.0
  n[2] = 0.0
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  FEMSei.Jacobi!(DF,detDF,pinvDF,X,Grids.Quad(),ksi1,ksi2,P1,P2,P3,P4,Rad)

  eInv = 1 / norm(pinvDF * n) 
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[2,1] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[2,2] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[2,3] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi4x(ksi1,ksi2) 
  grad[2] = gradphi4y(ksi1,ksi2) 
  ILoc[2,4] = eInv * (pinvDF * n)' * (pinvDF * grad)
  ksi1 = -1.0
  ksi2 = 1.0
  n[1] = 0.0
  n[2] = 1.0
  FEMSei.Jacobi!(DF,detDF,pinvDF,X,Grids.Quad(),ksi1,ksi2,P1,P2,P3,P4,Rad)
  eInv = 1 / norm(pinvDF * n) 
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[4,1] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[4,2] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[4,3] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi4x(ksi1,ksi2) 
  grad[2] = gradphi4y(ksi1,ksi2) 
  ILoc[4,4] = eInv * (pinvDF * n)' * (pinvDF * grad)
end

function LocalInterpolationGradTri!(ILoc,P1,P2,P3,n1,n2,Rad)

  grad = zeros(2)
  n = zeros(2)
  @polyvar x y ksi eta
  nu1 = -1.0*ksi - 1.0*eta + 1.0
  nu2 = 1.0*ksi + 0.0*eta + 0.0
  nu3 = 0.0*ksi + 1.0*eta + 0.0
  phi1 = subs(nu1, ksi => (x+1)/2, eta => (y+1)/2)
  phi2 = subs(nu2, ksi => (x+1)/2, eta => (y+1)/2)
  phi3 = subs(nu3, ksi => (x+1)/2, eta => (y+1)/2)

  gradphi1x = differentiate(phi1,x)
  gradphi1y = differentiate(phi1,y)
  gradphi2x = differentiate(phi2,x)
  gradphi2y = differentiate(phi2,y)
  gradphi3x = differentiate(phi3,x)
  gradphi3y = differentiate(phi3,y)
  @. ILoc = 0
  # Pointwise interpolation
  ILoc[1,1] = 1
  #First gradient

  ksi1 = 1.0
  ksi2 = -1.0
  n[1] = 0.0
  n[2] = 1.0
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  FEMSei.Jacobi!(DF,detDF,pinvDF,X,Grids.Tri(),ksi1,ksi2,P1,P2,P3,Rad)

  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[2,1] = n1' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[2,2] = n1' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[2,3] = n1' * (pinvDF * grad)
  ksi1 = -1.0
  ksi2 = 1.0
  n[1] = 1.0
  n[2] = 0.0
  FEMSei.Jacobi!(DF,detDF,pinvDF,X,Grids.Tri(),ksi1,ksi2,P1,P2,P3,Rad)
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[3,1] = -n2' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[3,2] = -n2' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[3,3] = -n2' * (pinvDF * grad)
end
function LocalInterpolationGradTriOld!(ILoc,P1,P2,P3,Rad)

  grad = zeros(2)
  n = zeros(2)
  @polyvar x y ksi eta
  nu1 = -1.0*ksi - 1.0*eta + 1.0
  nu2 = 1.0*ksi + 0.0*eta + 0.0
  nu3 = 0.0*ksi + 1.0*eta + 0.0
  phi1 = subs(nu1, ksi => (x+1)/2, eta => (y+1)/2)
  phi2 = subs(nu2, ksi => (x+1)/2, eta => (y+1)/2)
  phi3 = subs(nu3, ksi => (x+1)/2, eta => (y+1)/2)

  gradphi1x = differentiate(phi1,x)
  gradphi1y = differentiate(phi1,y)
  gradphi2x = differentiate(phi2,x)
  gradphi2y = differentiate(phi2,y)
  gradphi3x = differentiate(phi3,x)
  gradphi3y = differentiate(phi3,y)
  @. ILoc = 0
  # Pointwise interpolation
  ILoc[1,1] = 1
  #First gradient

  ksi1 = 1.0
  ksi2 = -1.0
  n[1] = 1.0
  n[2] = 0.0
  DF = zeros(3,2)
  detDF = zeros(1)
  pinvDF = zeros(3,2)
  X = zeros(3)
  FEMSei.Jacobi!(DF,detDF,pinvDF,X,Grids.Tri(),ksi1,ksi2,P1,P2,P3,Rad)

  eInv = 1 / norm(pinvDF * n) 
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[2,1] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[2,2] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[2,3] = eInv * (pinvDF * n)' * (pinvDF * grad)
  ksi1 = -1.0
  ksi2 = 1.0
  n[1] = 0.0
  n[2] = 1.0
  FEMSei.Jacobi!(DF,detDF,pinvDF,X,Grids.Tri(),ksi1,ksi2,P1,P2,P3,Rad)
  eInv = 1 / norm(pinvDF * n) 
  grad[1] = gradphi1x(ksi1,ksi2) 
  grad[2] = gradphi1y(ksi1,ksi2) 
  ILoc[3,1] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi2x(ksi1,ksi2) 
  grad[2] = gradphi2y(ksi1,ksi2) 
  ILoc[3,2] = eInv * (pinvDF * n)' * (pinvDF * grad)
  grad[1] = gradphi3x(ksi1,ksi2) 
  grad[2] = gradphi3y(ksi1,ksi2) 
  ILoc[3,3] = eInv * (pinvDF * n)' * (pinvDF * grad)
end
