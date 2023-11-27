function CubedGrid(backend,FT,n,OrientFace,Rad,nz,Topography) 

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  NumFaces=0
  NumGhostFaces=0
  Faces=Array{Face}(undef, 0)
  NumEdges=0
  Edges=Array{Edge}(undef, 0)
  NumNodes=0
  Nodes=Array{Node}(undef, 0)
  Form=""
  Type=""
  Dim=0
  Rad=0.0
  NumEdgesI=0
  NumEdgesB=0
  nBar3=zeros(0,0)
  nBar=zeros(0,0)
  colors=[[]]
  Spline_2d = Spline2D(zeros(0),zeros(0),zeros(0),0,0,0.0)
  BoundaryFaces = KernelAbstractions.zeros(backend,Int,0)
  InteriorFaces = KernelAbstractions.zeros(backend,Int,0)
return GridStruct{FT,
                     typeof(z),
                     typeof(BoundaryFaces)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumGhostFaces,
    Faces,
    NumEdges,
    Edges,
    NumNodes,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    NumEdgesI,
    NumEdgesB,
    nBar3,
    nBar,
    Topography,
    colors,
    Spline_2d,
    BoundaryFaces,
    InteriorFaces,
    )
end

