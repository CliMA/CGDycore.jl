struct SphericalAngleCart <: Metric  end

#spherical_angular_cart(a, b) = 1.0 - cos(atan(norm(cross(a,b))/dot(a,b)))
spherical_angular_cart(a, b) = acos(max(min(dot(a,b)/(norm(a)*norm(b)),1),-1))
(::SphericalAngleCart)(a, b) = spherical_angular_cart(a, b)

function interpolate(SrcGrid,DestGrid)
  Points1 = zeros(3,SrcGrid.NumNodes)
  for i = 1 : SrcGrid.NumNodes
    Points1[1,i] = SrcGrid.Nodes[i].P.x   
    Points1[2,i] = SrcGrid.Nodes[i].P.y   
    Points1[3,i] = SrcGrid.Nodes[i].P.z   
  end  

  tree1 = BallTree(Points1, SphericalAngleCart(), leafsize=10)

  Val = Array{Float64,1}(undef,0)
  RowInd = Array{Int,1}(undef,0)
  ColInd = Array{Int,1}(undef,0)
  for iFD = 1 : DestGrid.NumFaces
    Mid = SVector{3}(DestGrid.Faces[iFD].Mid.x,DestGrid.Faces[iFD].Mid.y,DestGrid.Faces[iFD].Mid.z)
    DestPolygon  = Grids.Polygon(DestGrid.Faces[iFD],DestGrid.Nodes) 
    DestTriangles = Grids.triangulate(DestPolygon)
    areaDest = Grids.area(DestTriangles)
    r = 0.0
    for N in DestGrid.Faces[iFD].N
      P = SVector{3}(DestGrid.Nodes[N].P.x,DestGrid.Nodes[N].P.y,DestGrid.Nodes[N].P.z)
      r = max(r, SphericalAngleCart()(Mid,P))
    end  
    idxsLnn, dist = nn(tree1, Mid)
#   r = max(r, 1.5*dist + 10.0 * EPS)
    r = max(r, 1.9999*dist + 10.0 * EPS)
    idxsL = inrange(tree1, Mid, r, true)
    FaceSource = zeros(Int,0)
    for N in idxsL
      for F in SrcGrid.Nodes[N].F
        push!(FaceSource,F)
      end
    end
    FaceSource = unique(FaceSource)
    for iFS in FaceSource
      areaLoc = Grids.intersect(DestGrid.Faces[iFD],DestGrid,SrcGrid.Faces[iFS],SrcGrid,2.0*EPS)
      push!(Val,areaLoc/areaDest)
      push!(RowInd,iFD)
      push!(ColInd,iFS)
    end  
  end
  Inter = sparse(RowInd,ColInd,Val)
end
