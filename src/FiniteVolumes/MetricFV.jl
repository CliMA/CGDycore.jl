function MetricFiniteVolume(backend,FT,Grid)
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes
  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  Rad = Grid.Rad
  PrimalVolume = KernelAbstractions.zeros(backend,FT,NumFaces)
  DualVolume = KernelAbstractions.zeros(backend,FT,NumNodes)
  DualVolumeCPU = KernelAbstractions.zeros(backend,FT,NumNodes)
  DualEdgeVolume = KernelAbstractions.zeros(backend,FT,2,2,NumEdges)
  PrimalEdge = KernelAbstractions.zeros(backend,FT,NumEdges)
  DualEdge = KernelAbstractions.zeros(backend,FT,NumEdges)

  nCPU = zeros(FT,3,NumEdges)
  n = KernelAbstractions.zeros(backend,FT,3,NumEdges)
  x = zeros(FT,3)
  k = zeros(FT,3)
  t = zeros(FT,3)
  for iE = 1 : Grid.NumEdges
    x[1] = Grid.Edges[iE].Mid.x
    x[2] = Grid.Edges[iE].Mid.y
    x[3] = Grid.Edges[iE].Mid.z
    k[1] = Grid.Edges[iE].Mid.x
    k[2] = Grid.Edges[iE].Mid.y
    k[3] = Grid.Edges[iE].Mid.z
    k = k / norm(k)
    t[1] = Grid.Edges[iE].t.x
    t[2] = Grid.Edges[iE].t.y
    t[3] = Grid.Edges[iE].t.z
    nCPU[:,iE] = cross(k,t) 
  end
  copyto!(n,nCPU)
  PrimalEdgeCPU = zeros(FT,NumEdges)
  DualEdgeCPU = zeros(FT,NumEdges)
  for iE = 1 : NumEdges
    PrimalEdgeCPU[iE] = Grids.SizeGreatCircle(Edges[iE],Nodes) * Grid.Rad
    iF1 = Edges[iE].F[1]
    iF2 = Edges[iE].F[2]
    Mid1 = Faces[iF1].Mid
    Mid2 = Faces[iF2].Mid
    DualEdgeCPU[iE] = Grids.SizeGreatCircle(Mid1,Mid2) * Grid.Rad
  end
  copyto!(PrimalEdge,PrimalEdgeCPU)

  PrimalVolumeCPU = zeros(FT,NumFaces)
  for iF = 1 : NumFaces
    PrimalVolumeCPU[iF] = Grids.AreaFace(Faces[iF],Nodes) * Grid.Rad^2
  end
  copyto!(PrimalVolume,PrimalVolumeCPU)

  PrimalCircum = Array{Grids.Point}(undef,NumFaces)
  for iF = 1 : NumFaces
    if length(Faces[iF].N) == 3  
      P1 = Nodes[Faces[iF].N[1]].P
      P2 = Nodes[Faces[iF].N[2]].P
      P3 = Nodes[Faces[iF].N[3]].P
      if Faces[iF].Orientation == 1
        PrimalCircum[iF] = Grids.CircumCenter(P1,P2,P3)
      else
        PrimalCircum[iF] = Grids.CircumCenter(P1,P3,P2)
      end
    else
      PrimalCircum[iF] = Faces[iF].Mid
    end  
  end
  DualEdgeVolumeCPU = zeros(FT,2,2,NumEdges)
  for iE = 1 : NumEdges
    P1 = Nodes[Edges[iE].N[1]].P
    P2 = Nodes[Edges[iE].N[2]].P
    PC1 = Faces[Edges[iE].F[1]].Mid
    PC2 = Faces[Edges[iE].F[2]].Mid
    PEM = Edges[iE].Mid
    DualEdgeVolumeCPU[1,1,iE] = Grids.AreaSphericalTriangle(P1,PEM,PC1) * Grid.Rad^2
    DualEdgeVolumeCPU[2,1,iE] = Grids.AreaSphericalTriangle(P2,PEM,PC1) * Grid.Rad^2
    DualEdgeVolumeCPU[1,2,iE] = Grids.AreaSphericalTriangle(P1,PEM,PC2) * Grid.Rad^2
    DualEdgeVolumeCPU[2,2,iE] = Grids.AreaSphericalTriangle(P2,PEM,PC2) * Grid.Rad^2
  end
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]
      if i == 1
        im1 = length(Grid.Faces[iF].N)
      else
        im1 = i - 1
      end
      iEm1 = Grid.Faces[iF].E[im1]
      iE = Grid.Faces[iF].E[i]
      DualVolumeCPU[iN] += Grids.AreaSphericalTriangle(Grid.Nodes[iN].P,Grid.Edges[iE].Mid,
        Grid.Faces[iF].Mid) * Rad^2
      DualVolumeCPU[iN] += Grids.AreaSphericalTriangle(Grid.Nodes[iN].P,Grid.Faces[iF].Mid,
        Grid.Edges[iEm1].Mid) * Rad^2
    end
  end

  copyto!(DualEdgeVolume,DualEdgeVolumeCPU)
  copyto!(DualVolume,DualVolumeCPU)
  copyto!(DualEdge,DualEdgeCPU)
  return MetricFiniteVolume{FT,
                    typeof(PrimalVolume),
                    typeof(n),
                    typeof(DualEdgeVolume)}(
    PrimalVolume,
    DualVolume,
    PrimalEdge,
    DualEdge,
    DualEdgeVolume,
    n,
    )
end

function KineticEnergy(K,uN,Metric,Grid)

  @. K = 0
  @inbounds for iE = 1 : Grid.NumEdges
    iF1  = Grid.Edges[iE].F[1]
    iF2  = Grid.Edges[iE].F[2]
    KLoc  = uN[iE] * uN[iE] 
    K[iF1] += KLoc * Metric.DualEdgeVolume[1,iE]
    K[iF2] += KLoc * Metric.DualEdgeVolume[2,iE]
  end    
  @. K /= Metric.PrimalVolume
end

function KineticEnergy(K,uN,uT,Metric,Grid)

  @. K = 0
  @inbounds for iE = 1 : Grid.NumEdges
    iF1  = Grid.Edges[iE].F[1]
    iF2  = Grid.Edges[iE].F[2]
    KLoc  = 0.5 *(uN[iE] * uN[iE] + uT[iE] * uT[iE])
    K[iF1] += KLoc * Metric.DualEdgeVolume[1,iE]
    K[iF2] += KLoc * Metric.DualEdgeVolume[2,iE]
  end
  @. K /= Metric.PrimalVolume
end

