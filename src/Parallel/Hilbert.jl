  #! Type for one entry in the 3D Hilbert curve description table.
  #!
  #! Values of the $mv item are:
  #! - 1/-1 means x+1/x-1
  #! - 2/-2 means y+1/y-1
  #! - 3/-3 means z+1/z-1
  mutable struct hilbert_table
    mv::Array{Int,1}
    next::Array{Int,1}
    num2ref::Array{Int,1}
    ref2num::Array{Int,1}
  end  

  function hilbert_table()
    mv = zeros(Int,0)
    next = zeros(Int,0)
    num2ref = zeros(Int,0)
    ref2num = zeros(Int,0)
    return hilbert_table(
      mv,
      next,
      num2ref,
      ref2num,
    )
  end  
  function fill_hilbert_table(mv::Array{Int,1},next::Array{Int,1},num2ref::Array{Int,1},ref2num::Array{Int,1})
    @. mv = mv 
    @. next = next 
    @. num2ref = num2ref
    @. ref2num = ref2num 
    return hilbert_table(
      mv,
      next,
      num2ref,
      ref2num,
    )
  end  
  function InitHilbert()
    tab = map(UInt64(1):UInt64(24)) do i
      hilbert_table()
    end  
    tab[ 1] = fill_hilbert_table([ 1, 3,-1, 2, 1,-3,-1],[ 1, 2, 2,20,20, 4, 4, 5],[0,1,5,4,6,7,3,2],[0,1,7,6,3,2,4,5])  
    tab[ 2] = fill_hilbert_table([ 3, 2,-3, 1, 3,-2,-3],[ 6, 7, 7,19,19, 9, 9,10],[0,4,6,2,3,7,5,1],[0,7,3,4,1,6,2,5])  
    tab[ 3] = fill_hilbert_table([ 1, 2,-1, 3, 1,-2,-1],[11, 0, 0,23,23,14,14,15],[0,1,3,2,6,7,5,4],[0,1,3,2,7,6,4,5])  
    tab[ 4] = fill_hilbert_table([-1,-3, 1, 2,-1, 3, 1],[16,17,17, 7, 7,13,13, 8],[5,4,0,1,3,2,6,7],[2,3,5,4,1,0,6,7])  
    tab[ 5] = fill_hilbert_table([ 1,-2,-1,-3, 1, 2,-1],[19,14,14,10,10, 0, 0,22],[6,7,5,4,0,1,3,2],[4,5,7,6,3,2,0,1])  
    tab[ 6] = fill_hilbert_table([ 3,-2,-3,-1, 3, 2,-3],[23, 9, 9,15,15, 7, 7,21],[3,7,5,1,0,4,6,2],[4,3,7,0,5,2,6,1])  
    tab[ 7] = fill_hilbert_table([ 2, 1,-2, 3, 2,-1,-2],[ 0,11,11,13,13,15,15,14],[0,2,3,1,5,7,6,4],[0,3,1,2,7,4,6,5])  
    tab[ 8] = fill_hilbert_table([ 3, 1,-3, 2, 3,-1,-3],[ 2, 1, 1, 3, 3, 5, 5, 4],[0,4,5,1,3,7,6,2],[0,3,7,4,1,2,6,5])  
    tab[ 9] = fill_hilbert_table([-3,-2, 3, 1,-3, 2, 3],[21,18,18,11,11,20,20,23],[6,2,0,4,5,1,3,7],[2,5,1,6,3,4,0,7])  
    tab[10] = fill_hilbert_table([ 3,-1,-3,-2, 3, 1,-3],[13, 5, 5,14,14, 1, 1,17],[3,7,6,2,0,4,5,1],[4,7,3,0,5,6,2,1])  
    tab[11] = fill_hilbert_table([-3,-2, 3, 1,-3, 2, 3],[ 3,15,15, 4, 4,11,11,12],[5,7,6,4,0,2,3,1],[4,7,5,6,3,0,2,1])  
    tab[12] = fill_hilbert_table([ 3,-1,-3,-2, 3, 1,-3],[ 7, 6, 6, 8, 8,10,10, 9],[0,2,6,4,5,7,3,1],[0,7,1,6,3,4,2,5])  
    tab[13] = fill_hilbert_table([ 2,-1,-2,-3, 2, 1,-2],[ 5,13,13,18,18,17,17, 1],[3,2,6,7,5,4,0,1],[6,7,1,0,5,4,2,3])  
    tab[14] = fill_hilbert_table([ 2, 3,-2, 1, 2,-3,-2],[22,12,12, 6, 6, 3, 3,19],[3,2,0,1,5,4,6,7],[2,3,1,0,5,4,6,7])  
    tab[15] = fill_hilbert_table([-1, 3, 1,-2,-1,-3, 1],[ 8, 4, 4, 9, 9, 2, 2,16],[6,7,3,2,0,1,5,4],[4,5,3,2,7,6,0,1])  
    tab[16] = fill_hilbert_table([-1,-2, 1, 3,-1, 2, 1],[20,10,10, 5, 5, 6, 6,18],[5,7,3,1,0,2,6,4],[4,3,5,2,7,0,6,1])  
    tab[17] = fill_hilbert_table([ 1,-3,-1,-2, 1, 3,-1],[10,20,20,22,22,18,18, 6],[5,1,3,7,6,2,0,4],[6,1,5,2,7,0,4,3])  
    tab[18] = fill_hilbert_table([ 2,-3,-2,-1, 2, 3,-2],[15, 3, 3,21,21,12,12,11],[5,4,6,7,3,2,0,1],[6,7,5,4,1,0,2,3])  
    tab[19] = fill_hilbert_table([-3, 2, 3,-1,-3,-2, 3],[ 4, 8, 8,12,12,16,16, 2],[6,2,3,7,5,1,0,4],[6,5,1,2,7,4,0,3])  
    tab[20] = fill_hilbert_table([-1, 2, 1,-3,-1,-2, 1],[18,21,21, 1, 1,23,23,20],[6,4,0,2,3,1,5,7],[2,5,3,4,1,6,0,7])  
    tab[21] = fill_hilbert_table([-3, 1, 3,-2,-3,-1, 3],[17,16,16, 0, 0, 8, 8,13],[5,1,0,4,6,2,3,7],[2,1,5,6,3,0,4,7])  
    tab[22] = fill_hilbert_table([-2, 1, 2,-3,-2,-1, 2],[14,19,19,17,17,22,22, 0],[6,4,5,7,3,1,0,2],[6,5,7,4,1,2,0,3])  
    tab[23] = fill_hilbert_table([-2, 3, 2,-1,-2,-3, 2],[ 9,23,23,16,16,21,21, 7],[3,1,5,7,6,4,0,2],[6,1,7,0,5,2,4,3])  
    tab[24] = fill_hilbert_table([-2,-1, 2, 3,-2, 1, 2],[12,22,22, 2, 2,19,19, 3],[3,1,0,2,6,4,5,7],[2,1,3,0,5,6,4,7])  
    return tab
  end
    
  

function HilbertOrderXY(Point,P0,P1,lev::Int,tab)

  HilbertOrderXY = 1
  s = 0
  dim=2^lev
  # length of curve in base octant
  len=dim^3
  PMin = deepcopy(P0)
  PMax = deepcopy(P1)
  PCenter = SVector{3}(0.5*(PMin+PMax))

  for l = lev:-1:1
    # dimension and length of refined octants
    dim=dim/2
    len=len/8
    
    # determine spatial octant and reduce spatial coordinates
    # so that they are relative to refined octant
    # reference octant:
    # z=0 | z=1  
    # 2 3 | 6 7
    # 0 1 | 4 5
    rq = UInt64(0)
    OneInt8 = UInt64(1)
    TwoInt8 = UInt64(2)
    if Point[1] >= PCenter[1] 
      rq = OneInt8 | rq
      PMin[1] = PCenter[1]
    else  
      PMax[1] = PCenter[1]
    end
    if Point[2] >= PCenter[2]
      rq = TwoInt8 | rq
      PMin[2] = PCenter[2]
    else  
      PMax[2] = PCenter[2]
    end
    
    # get sequence number of octant q in the curve primitive

    q = tab[s+1].ref2num[rq+1]
    
    # add offset to index    
    HilbertOrderXY = HilbertOrderXY + q * len
    
    # get next state
    s = tab[s+1].next[q+1]
    PCenter = SVector{3}(0.5 * (PMin + PMax))
  end
  return HilbertOrderXY
end

function HilbertFaceSphere!(Grid,P0Sph,P1Sph)
  hilbert_table = InitHilbert()
  lev = 10
  FaceOrder = zeros(Int,Grid.NumFaces)
  PointSph = zeros(3)
  for iF = 1 : Grid.NumFaces
    Point = Grid.Faces[iF].Mid
    (PointSph[1],PointSph[2],PointSph[3])=CGDycore.cart2sphere(Point.x,Point.y,Point.z)
    FaceOrder[iF] = HilbertOrderXY(PointSph,P0Sph,P1Sph,lev,hilbert_table)
  end
  p = sortperm(FaceOrder)
  permute!(Grid.Faces,p)
  for iF = 1 : Grid.NumFaces
    FaceOrder[Grid.Faces[iF].F] = iF
    Grid.Faces[iF].F = iF
  end
  for iE = 1 : Grid.NumEdges
     if Grid.Edges[iE].F[1] > 0
       Grid.Edges[iE].F[1] = FaceOrder[Grid.Edges[iE].F[1]]
     end
     if Grid.Edges[iE].F[2] > 0
       Grid.Edges[iE].F[2] = FaceOrder[Grid.Edges[iE].F[2]]
     end
  end
  for iN = 1 : Grid.NumNodes
    for i in eachindex(Grid.Nodes[iN].F)
      Grid.Nodes[iN].F[i] = FaceOrder[Grid.Nodes[iN].F[i]]
    end
  end
end  
