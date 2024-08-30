function SphericalGrid(backend,FT,nLon,nLat,LatB,OrientFace,Rad,nz)

  nBar=[ 0  1   0   1
        -1  0  -1   0];
  Dim=3;
  Type=Quad()
  Rad=Rad;
  Form="Sphere";
  
  NumNodes = nLon * (nLat + 1)
  NumNodesB = nLon * 2
  Nodes = map(1:NumNodes) do i
    Node()
  end
  dLon = 2 * pi / nLon
  dLat = 2 * LatB /nLat
  P=zeros(Float64,3,nLon+1,nLat+1)
  Lat = -LatB
  @inbounds for iLat=1:nLat+1
    Lon = 0.0
    @inbounds for iLon=1:nLon+1
      P[:,iLon,iLat] .= sphere2cart(Lon,Lat,Rad)
      Lon = Lon + dLon
    end
    Lat = Lat + dLat
  end
  NodeNumber = 1
  Lon = 0.0
  Lat = -LatB  
  for i = 1 : nLat + 1 
    Lon = 0.0
    for j = 1 : nLon
      PP = sphere2cart(Lon,Lat,Rad)
      if i == 1 || i == nLat + 1
        Nodes[NodeNumber] = Node(Point(PP),NodeNumber,'B')
      else  
        Nodes[NodeNumber] = Node(Point(PP),NodeNumber,' ')
      end  
      NodeNumber += 1
      Lon += dLon
    end
    Lat += dLat
  end  


  NumEdges=nLon*nLat+nLon*(nLat+1)
  NumEdgesLon=nLon*(nLat+1)

  Edges= map(1:NumEdges) do i
    Edge([1,2],Nodes,0,0,"",0)
  end

  EdgeNumber=1
  EdgeNumberLon=1
  EdgeNumberLat=1
  BC=""
  N1=1
  N2=nLon+1
  @inbounds for iLat=1:nLat
    @inbounds for iLon=1:nLon
      Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"",EdgeNumberLat)
      EdgeNumber=EdgeNumber+1
      EdgeNumberLat=EdgeNumberLat+1
      N1=N1+1
      N2=N2+1
    end
  end
  N1=1
  N2=2
  @inbounds for iLat=1:nLat+1
    if iLat == 1 || iLat == nLat+1
      TypeE = "B"  
    else
      TypeE = ""
    end  
    @inbounds for iLon=1:nLon
      if iLon == nLon 
        Edges[EdgeNumber] = Edge([N1,1+(iLat-1)*nLon],Nodes,EdgeNumber,EdgeNumber,TypeE,EdgeNumberLon)
        EdgeNumber = EdgeNumber + 1
        EdgeNumberLon = EdgeNumberLon + 1
        N1 += 1
        N2 += 1
      else
        Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,TypeE,EdgeNumberLon)
        EdgeNumber=EdgeNumber+1
        EdgeNumberLon=EdgeNumberLon+1
        N1=N1+1
        N2=N2+1
        if iLon==nLon
          N1=N1+1
          N2=N2+1
        end
      end
    end
  end

  NumFaces=nLon*nLat
  Faces=map(1:NumFaces) do i
     Face()
  end
  E1=nLon*nLat+1
  E3=E1+nLon
  E2=2
  E4=1
  FaceNumber=1
  TypeF="o"
  @inbounds for iLat=1:nLat
    @inbounds for iLon=1:nLon
      if iLon == nLon 
        (Faces[FaceNumber],Edges)=Face([E1,1+(iLat-1)*nLon,E3,E4],Nodes,Edges,FaceNumber,TypeF,OrientFace,
          P=[P[:,iLon,iLat] P[:,iLon+1,iLat] P[:,iLon+1,iLat+1] P[:,iLon,iLat+1]])
        FaceNumber=FaceNumber+1
        E1=E1+1
        E3=E3+1
      else  
        (Faces[FaceNumber],Edges)=Face([E1,E2,E3,E4],Nodes,Edges,FaceNumber,TypeF,OrientFace,
          P=[P[:,iLon,iLat] P[:,iLon+1,iLat] P[:,iLon+1,iLat+1] P[:,iLon,iLat+1]])
        FaceNumber=FaceNumber+1
        E1=E1+1
        E2=E2+1
        E3=E3+1
        E4=E4+1
      end
    end
    E2=E2+1
    E4=E4+1
  end

  Orientation!(Edges,Faces);
  Renumbering!(Edges,Faces);
  FacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  nBar3 = zeros(0,0)
  nBar = zeros(0,0)
  NumBoundaryFaces = 0
  NumFacesB = 0
  NumFacesG = 0
  NumEdgesB = 0
  NumEdgesG = 0
  NumNodesG = 0
  AdaptGrid = ""
  return GridStruct{FT,
                    typeof(z)}(
    nz,
    zP,
    z,
    dzeta,
    H,
    NumFaces,
    NumFacesB,
    NumFacesG,
    Faces,
    NumEdges,
    NumEdgesB,
    NumEdgesG,
    Edges,
    NumNodes,
    NumNodesB,
    NumNodesG,
    Nodes,
    Form,
    Type,
    Dim,
    Rad,
    nBar3,
    nBar,
    AdaptGrid,
    )

end
