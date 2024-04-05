function MetricFiniteVolume(backend,FT,Grid,type::Grids.Tri)
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes
  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  PrimalVolume = KernelAbstractions.zeros(backend,FT,NumFaces)
  DualVolume = KernelAbstractions.zeros(backend,FT,0)
  DualEdgeVolume = KernelAbstractions.zeros(backend,FT,NumEdges)
  PrimalEdge = KernelAbstractions.zeros(backend,FT,NumEdges)
  DualEdge = KernelAbstractions.zeros(backend,FT,NumEdges)

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
    P1 = Nodes[Faces[iF].N[1]].P
    P2 = Nodes[Faces[iF].N[2]].P
    P3 = Nodes[Faces[iF].N[3]].P
    if Faces[iF].Orientation == 1
      PrimalCircum[iF] = Grids.CircumCenter(P1,P2,P3)
    else
      PrimalCircum[iF] = Grids.CircumCenter(P1,P3,P2)
    end
  end
  DualEdgeVolumeCPU = zeros(FT,NumEdges)
  for iE = 1 : NumEdges
    P1 = Nodes[Edges[iE].N[1]].P
    P2 = Nodes[Edges[iE].N[2]].P
    PC1 = PrimalCircum[Edges[iE].F[1]]
    PC2 = PrimalCircum[Edges[iE].F[2]]
    PC1 = Faces[Edges[iE].F[1]].Mid
    PC2 = Faces[Edges[iE].F[2]].Mid
    DualEdgeVolumeCPU[iE] = (Grids.AreaSphericalTriangle(P1,P2,PC1) +
      Grids.AreaSphericalTriangle(P1,P2,PC2)) * Grid.Rad^2
  end
  copyto!(DualEdgeVolume,DualEdgeVolumeCPU)
  copyto!(DualEdge,DualEdgeCPU)
  return MetricFiniteVolume{FT,
                    typeof(PrimalVolume)}(
    PrimalVolume,
    DualVolume,
    PrimalEdge,
    DualEdge,
    DualEdgeVolume,
    )
end
