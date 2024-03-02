function CartGrid(backend,FT,nx::Int,ny::Int,lx::Float64,ly::Float64,x0::Float64,y0::Float64,OrientFace,Boundary,nz)
  nBar=[ 0  1   0   1
             -1  0  -1   0]
  nBar3=[ 0  1   0   1
             -1  0  -1   0
             0   0  0    0]
  Dim=3
  Type = Quad()
  Form="Planar"
  Pert=0.0
  PertX=0.2
  PertY=0.2
  if Boundary.WE == "Period" && Boundary.SN == "Period"
    NumNodes=nx*ny
  elseif Boundary.WE == "Period"
    NumNodes=nx*(ny+1)
  elseif Boundary.SN == "Period"
    NumNodes=(nx+1)*ny
  else
    NumNodes=(nx+1)*(ny+1)
  end
  Nodes = map(1:NumNodes) do i
    Node()
  end
  NodeNumber=1
  dx=lx/nx
  dy=ly/ny
  y=y0
  P=zeros(Float64,3,nx+1,ny+1)
  @inbounds for iy=1:ny+1
    eta=(iy-1)/ny
    x=x0
    @inbounds for ix=1:nx+1
      P[1,ix,iy]=x
      P[2,ix,iy]=y
      P[3,ix,iy]=0
      x=x+dx
    end
    y=y+dy
  end

  y=y0
  @inbounds for iy=1:ny+1
    x=x0
    if iy==ny+1 && Boundary.SN == "Period"
    else
      @inbounds for ix=1:nx+1
        if ix==nx+1 && Boundary.WE == "Period"
        else
          Nodes[NodeNumber]=Node(Point([x,y,0.0]),NodeNumber,' ')
          NodeNumber=NodeNumber+1
        end
        x=x+dx
      end
    end
    y=y+dy
  end

  if Boundary.WE == "Period"  && Boundary.SN == "Period"
    NumEdges=2*nx*ny
    NumEdgesX=nx*ny
  elseif Boundary.WE == "Period"
    NumEdges=nx*ny+nx*(ny+1)
    NumEdgesX=nx*(ny+1)
  elseif Boundary.SN == "Period"
    NumEdges=(nx+1)*ny+nx*ny
    NumEdgesX=nx*ny
  else
    NumEdges=(nx+1)*ny+nx*(ny+1)
    NumEdgesX=(nx+1)*ny
  end

  Edges= map(1:NumEdges) do i
    Edge([1,2],Nodes,0,0,"",0)
  end

  EdgeNumber=1
  EdgeNumberX=1
  EdgeNumberY=1
  BC=""
  if Boundary.WE == "Period" && Boundary.SN == "Period"
    N1=1
    N2=nx+1
  elseif Boundary.WE == "Period"
    N1=1
    N2=nx+1
  elseif Boundary.SN == "Period"
    N1=1
    N2=nx+2
  else
    N1=1
    N2=nx+2
  end
  @inbounds for iy=1:ny
    @inbounds for ix=1:nx+1
      if ix==nx+1 && Boundary.WE == "Period"
      else
        if iy==ny && Boundary.SN == "Period"
          Edges[EdgeNumber]=Edge([N1,1+(ix-1)],Nodes,EdgeNumber,EdgeNumber,"Y",EdgeNumberY)
          EdgeNumber=EdgeNumber+1
          EdgeNumberY=EdgeNumberY+1
          N1=N1+1
          N2=N2+1
        else
          Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,"Y",EdgeNumberY)
          EdgeNumber=EdgeNumber+1
          EdgeNumberY=EdgeNumberY+1
          N1=N1+1
          N2=N2+1
        end
      end
    end
  end
  N1=1
  N2=2
  TypeE = "X"
  @inbounds for iy=1:ny+1
    if iy==ny+1 && Boundary.SN == "Period"
    else
      @inbounds for ix=1:nx
        if ix==nx && Boundary.WE == "Period"
          Edges[EdgeNumber]=Edge([N1,1+(iy-1)*nx],Nodes,EdgeNumber,EdgeNumber,TypeE,EdgeNumberX)
          EdgeNumber=EdgeNumber+1
          EdgeNumberX=EdgeNumberX+1
          N1=N1+1
          N2=N2+1
        else
          Edges[EdgeNumber]=Edge([N1,N2],Nodes,EdgeNumber,EdgeNumber,TypeE,EdgeNumberX)
          EdgeNumber=EdgeNumber+1
          EdgeNumberX=EdgeNumberX+1
          N1=N1+1
          N2=N2+1
          if ix==nx
            N1=N1+1
            N2=N2+1
          end
        end
      end
    end
  end

  NumFaces=nx*ny
  Faces=map(1:NumFaces) do i
     Face()
  end
  if Boundary.WE == "Period" && Boundary.SN == "Period"
    E1=nx*ny+1
    E3=nx*ny+1+nx
  else
    E1=(nx+1)*ny+1
    E3=(nx+1)*ny+1+nx
  end
  E2=2
  E4=1
  FaceNumber=1
  TypeF="o"
  @inbounds for iy=1:ny
    if iy==ny && Boundary.SN == "Period"
      @inbounds for ix=1:nx
        if ix==nx && Boundary.WE == "Period"
          (Faces[FaceNumber],Edges)=Face([E1,1+(iy-1)*nx,NumEdgesX+1+(ix-1),E4],Nodes,Edges,FaceNumber,TypeF,OrientFace,
            P=[P[:,ix,iy] P[:,ix+1,iy] P[:,ix+1,iy+1] P[:,ix,iy+1]])
          FaceNumber=FaceNumber+1
          E1=E1+1
          E4=E4+1
        else
          (Faces[FaceNumber],Edges)=Face([E1,E2,NumEdgesX+1+(ix-1),E4],Nodes,Edges,FaceNumber,TypeF,OrientFace,
            P=[P[:,ix,iy] P[:,ix+1,iy] P[:,ix+1,iy+1] P[:,ix,iy+1]])
          FaceNumber=FaceNumber+1
          E1=E1+1
          E2=E2+1
          E4=E4+1
        end
      end
    else
      @inbounds for ix=1:nx
        if ix==nx && Boundary.WE == "Period"
          (Faces[FaceNumber],Edges)=Face([E1,1+(iy-1)*nx,E3,E4],Nodes,Edges,FaceNumber,TypeF,OrientFace,
            P=[P[:,ix,iy] P[:,ix+1,iy] P[:,ix+1,iy+1] P[:,ix,iy+1]])
          FaceNumber=FaceNumber+1
          E1=E1+1
          E3=E3+1
        else
          (Faces[FaceNumber],Edges)=Face([E1,E2,E3,E4],Nodes,Edges,FaceNumber,TypeF,OrientFace,
            P=[P[:,ix,iy] P[:,ix+1,iy] P[:,ix+1,iy+1] P[:,ix,iy+1]])
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
  end


  Orientation!(Edges,Faces);
  Renumbering!(Edges,Faces);
  FacesInNodes!(Nodes,Faces)

  zP=zeros(nz)
  z=KernelAbstractions.zeros(backend,FT,nz+1)
  dzeta=zeros(nz)
  H=0.0
  colors=[[]]
  NumGhostFaces = 0
  nBar3 = zeros(0,0)
  nBar = zeros(0,0)
  NumBoundaryFaces = 0
  NumEdgesI = 0
  NumEdgesB = 0
  AdaptGrid = ""
  return GridStruct{FT,
                    typeof(z)}(
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
    colors,
    NumBoundaryFaces,
    AdaptGrid,
    )

end
