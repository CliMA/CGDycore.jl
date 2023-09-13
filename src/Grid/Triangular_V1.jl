module Triangular

struct Point{FT}
  x::FT
  y::FT
  z::FT
end  

struct NodeTri_T{FT}
  P::Point{FT}
  Number::Int
  Edge::Array{Int,1}
  Next::Union{NodeTri_T{FT}, Nothing}
end

mutable struct LinkedNodeTri_t{T}
    head::Union{NodeTri_{T}, Nothing}
end

LinkedList(T::Type) = LinkedNodeTri_T{T}(nothing)

Base.iterate(ll::LinkedNodeTri_T) = ll.head === nothing ? nothing : (ll.head.value, ll.head)
Base.iterate(ll::LinkedNodeTri_T{T}, state::NodeTri_T{T}) where T =
    state.next === nothing ? nothing : (state.next.value, state.next)

function Base.getindex(ll::LinkedNodeTri_T, idx::Integer)
    idx < 1 && throw(BoundsError("$idx is less than 1"))
    for v in ll
        idx -= 1
        idx == 0 && return v
    end
    throw(BoundsError("index beyond end of linked list"))
end

function Base.pushfirst!(ll::LinkedNodeTri_T{T}, items::T...) where T
    for item in reverse(items)
        ll.head = NodeTri_TNode{T}(item, ll.head)
    end
    ll
end

Base.show(io::IO, ll::LinkedNodeTri_T{T}) where T =
    print(io, "LinkedNodeTri_T{$T}[" * join(ll, ", ") * "]")

Base.eltype(ll::LinkedNodeTri_T{T}) where T = T

Base.length(ll::LinkedNodeTri_T) = count(v -> true, ll)

Base.firstindex(ll::LinkedNodeTri_T) = 1
Base.lastindex(ll::LinkedNodeTri_T) = length(ll)

end

#=
MODULE Triangular_Mod

 USE Geometry_Mod
 USE Control_Mod

 IMPLICIT NONE

 TYPE NodeTri_T
   TYPE(Point_T) :: P
   INTEGER :: Number
   TYPE(NodeTri_T), POINTER :: Next=>NULL()
   INTEGER, POINTER :: Edge(:)
 END TYPE NodeTri_T

 TYPE NodeTriVec_T
   TYPE(NodeTri_T), POINTER :: NodeTri=>NULL()
 END TYPE NodeTriVec_T

 TYPE EdgeTri_T
   INTEGER :: Number
   TYPE(NodeTri_T), POINTER :: Node1,Node2
   INTEGER, POINTER :: Cell(:)
   TYPE(EdgeTri_T), POINTER :: Next=>NULL()
 END TYPE EdgeTri_T

 TYPE EdgeTriVec_T
   TYPE(EdgeTri_T), POINTER :: EdgeTri=>NULL()
 END TYPE EdgeTriVec_T

 TYPE CellTri_T
   INTEGER :: Number
   TYPE(EdgeTri_T), POINTER :: Edge1,Edge2,Edge3
   INTEGER :: OrientE1=1,OrientE2=1,OrientE3=1
   TYPE(CellTri_T), POINTER :: Next=>NULL()
 END TYPE CellTri_T

 TYPE CellTriVec_T
   TYPE(CellTri_T), POINTER :: CellTri=>NULL()
 END TYPE CellTriVec_T

 TYPE TriangularGrid_T
   TYPE(NodeTri_T), POINTER :: StartNode=>NULL()
   TYPE(EdgeTri_T), POINTER :: StartEdge=>NULL()
   TYPE(CellTri_T), POINTER :: StartCell=>NULL()
   TYPE(NodeTri_T), POINTER :: CurrentNode=>NULL()
   TYPE(EdgeTri_T), POINTER :: CurrentEdge=>NULL()
   TYPE(CellTri_T), POINTER :: CurrentCell=>NULL()
   TYPE(CellTri_T), POINTER :: LastCell=>NULL()
   TYPE(NodeTri_T), POINTER :: LastNode=>NULL()
   TYPE(EdgeTri_T), POINTER :: LastEdge=>NULL()
 END TYPE TriangularGrid_T

 INTERFACE OPERATOR(.INS.)
   MODULE PROCEDURE AppendNode,AppendEdge,AppendCell
 END INTERFACE
 INTERFACE CountNode
   MODULE PROCEDURE CountNode
 END INTERFACE  
 INTERFACE CountEdge
   MODULE PROCEDURE CountEdge
 END INTERFACE  
 INTERFACE CountCell
   MODULE PROCEDURE CountCell
 END INTERFACE  
 INTERFACE MidPoint
   MODULE PROCEDURE MidPoint
 END INTERFACE  
 INTERFACE SquareGrid
   MODULE PROCEDURE SquareGrid
 END INTERFACE  
 INTERFACE Deallocate
   MODULE PROCEDURE DeallocateTriangularGrid
 END INTERFACE

CONTAINS

SUBROUTINE SquareGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid

  TYPE(CellTri_T), POINTER :: CellTri1
  TYPE(CellTri_T), POINTER :: CellTri2
  TYPE(EdgeTri_T), POINTER :: EdgeTri1
  TYPE(EdgeTri_T), POINTER :: EdgeTri2
  TYPE(EdgeTri_T), POINTER :: EdgeTri3
  TYPE(EdgeTri_T), POINTER :: EdgeTri4
  TYPE(EdgeTri_T), POINTER :: EdgeTri5
  TYPE(NodeTri_T), POINTER :: NodeTri1
  TYPE(NodeTri_T), POINTER :: NodeTri2
  TYPE(NodeTri_T), POINTER :: NodeTri3
  TYPE(NodeTri_T), POINTER :: NodeTri4

  ALLOCATE(CellTri1)
  ALLOCATE(CellTri2)
  ALLOCATE(EdgeTri1)
  ALLOCATE(EdgeTri2)
  ALLOCATE(EdgeTri3)
  ALLOCATE(EdgeTri4)
  ALLOCATE(EdgeTri5)
  ALLOCATE(NodeTri1)
  ALLOCATE(NodeTri2)
  ALLOCATE(NodeTri3)
  ALLOCATE(NodeTri4)
  NodeTri1%P%x=0.0d0
  NodeTri1%P%y=0.0d0
  NodeTri1%P%z=0.0d0
  NodeTri2%P%x=1.0d0
  NodeTri2%P%y=0.0d0
  NodeTri2%P%z=0.0d0
  NodeTri3%P%x=1.0d0
  NodeTri3%P%y=1.0d0
  NodeTri3%P%z=0.0d0
  NodeTri4%P%x=0.0d0
  NodeTri4%P%y=1.0d0
  NodeTri4%P%z=0.0d0
  EdgeTri1%Node1=>NodeTri1
  EdgeTri1%Node2=>NodeTri2
  EdgeTri2%Node1=>NodeTri2
  EdgeTri2%Node2=>NodeTri3
  EdgeTri3%Node1=>NodeTri3
  EdgeTri3%Node2=>NodeTri4
  EdgeTri4%Node1=>NodeTri4
  EdgeTri4%Node2=>NodeTri1
  EdgeTri5%Node1=>NodeTri1
  EdgeTri5%Node2=>NodeTri3
  CellTri1%Edge1=>EdgeTri1
  CellTri1%OrientE1=1
  CellTri1%Edge2=>EdgeTri2
  CellTri1%OrientE2=1
  CellTri1%Edge3=>EdgeTri5
  CellTri1%OrientE3=-1
  CellTri2%Edge1=>EdgeTri3
  CellTri2%OrientE1=1
  CellTri2%Edge2=>EdgeTri4
  CellTri2%OrientE2=1
  CellTri2%Edge3=>EdgeTri5
  CellTri2%OrientE3=1
  CALL InsertCell(CellTri1,TriangularGrid)
  CALL InsertCell(CellTri2,TriangularGrid)
  CALL InsertEdge(EdgeTri1,TriangularGrid)
  CALL InsertEdge(EdgeTri2,TriangularGrid)
  CALL InsertEdge(EdgeTri3,TriangularGrid)
  CALL InsertEdge(EdgeTri4,TriangularGrid)
  CALL InsertEdge(EdgeTri5,TriangularGrid)
  CALL InsertNode(NodeTri1,TriangularGrid)
  CALL InsertNode(NodeTri2,TriangularGrid)
  CALL InsertNode(NodeTri3,TriangularGrid)
  CALL InsertNode(NodeTri4,TriangularGrid)
END SUBROUTINE SquareGrid

SUBROUTINE IcosahedronGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid

  INTEGER :: i
  REAL(8) :: lam,phi,Pi
  TYPE(NodeTri_T), POINTER :: NodeTop,NodeBottom
  TYPE(NodeTriVec_T), POINTER :: NodeLayerTop(:),NodeLayerBottom(:)
  TYPE(EdgeTriVec_T), POINTER :: EdgeTop(:),EdgeBottom(:)
  TYPE(EdgeTriVec_T), POINTER :: EdgeLayerTop(:),EdgeLayerBottom(:)
  TYPE(EdgeTriVec_T), POINTER :: EdgeMid(:)
  TYPE(CellTriVec_T), POINTER :: CellTop(:),CellBottom(:)
  TYPE(CellTriVec_T), POINTER :: CellMid(:)

  Pi=4.0d0*ATAN(1.0d0)

  ALLOCATE(NodeLayerTop(5),NodeLayerBottom(5))
  ALLOCATE(EdgeTop(5),EdgeBottom(5))
  ALLOCATE(EdgeLayerTop(5),EdgeLayerBottom(5))
  ALLOCATE(EdgeMid(10))
  ALLOCATE(CellTop(5),CellBottom(5))
  ALLOCATE(CellMid(10))

  ALLOCATE(NodeTop)
  NodeTop%P%x=0.0d0
  NodeTop%P%y=0.0d0
  NodeTop%P%z=RadOut
  CALL InsertNode(NodeTop,TriangularGrid)
  phi=ATAN(0.5d0)
  lam=0.0d0
  DO i=1,5
    ALLOCATE(NodeLayerTop(i)%NodeTri)
    NodeLayerTop(i)%NodeTri%P%x=RadOut*COS(lam)*COS(phi)
    NodeLayerTop(i)%NodeTri%P%y=RadOut*SIN(lam)*COS(phi)
    NodeLayerTop(i)%NodeTri%P%z=RadOut*SIN(phi)
    lam=lam+2.0d0*Pi/5.0d0
    CALL InsertNode(NodeLayerTop(i)%NodeTri,TriangularGrid)
  END DO  
  phi=-ATAN(0.5d0)
  lam=Pi/5.0d0
  DO i=1,5
    ALLOCATE(NodeLayerBottom(i)%NodeTri)
    NodeLayerBottom(i)%NodeTri%P%x=RadOut*COS(lam)*COS(phi)
    NodeLayerBottom(i)%NodeTri%P%y=RadOut*SIN(lam)*COS(phi)
    NodeLayerBottom(i)%NodeTri%P%z=RadOut*SIN(phi)
    lam=lam+2.0d0*Pi/5.0d0
    CALL InsertNode(NodeLayerBottom(i)%NodeTri,TriangularGrid)
  END DO  
  ALLOCATE(NodeBottom)
  NodeBottom%P%x=0.0d0
  NodeBottom%P%y=0.0d0
  NodeBottom%P%z=-RadOut
  CALL InsertNode(NodeBottom,TriangularGrid)
! Edges  
  DO i=1,5
    ALLOCATE(EdgeTop(i)%EdgeTri)
    EdgeTop(i)%EdgeTri%Node1=>NodeLayerTop(i)%NodeTri
    EdgeTop(i)%EdgeTri%Node2=>NodeTop
    CALL InsertEdge(EdgeTop(i)%EdgeTri,TriangularGrid)
  END DO  
  DO i=1,4
    ALLOCATE(EdgeLayerTop(i)%EdgeTri)
    EdgeLayerTop(i)%EdgeTri%Node1=>NodeLayerTop(i)%NodeTri
    EdgeLayerTop(i)%EdgeTri%Node2=>NodeLayerTop(i+1)%NodeTri
    CALL InsertEdge(EdgeLayerTop(i)%EdgeTri,TriangularGrid)
  END DO  
  ALLOCATE(EdgeLayerTop(5)%EdgeTri)
  EdgeLayerTop(5)%EdgeTri%Node1=>NodeLayerTop(5)%NodeTri
  EdgeLayerTop(5)%EdgeTri%Node2=>NodeLayerTop(1)%NodeTri
  CALL InsertEdge(EdgeLayerTop(5)%EdgeTri,TriangularGrid)
  DO i=1,10
    ALLOCATE(EdgeMid(i)%EdgeTri)
    CALL InsertEdge(EdgeMid(i)%EdgeTri,TriangularGrid)
  END DO  
  DO i=1,4
    EdgeMid(2*i-1)%EdgeTri%Node1=>NodeLayerBottom(i)%NodeTri
    EdgeMid(2*i-1)%EdgeTri%Node2=>NodeLayerTop(i)%NodeTri
    EdgeMid(2*i)%EdgeTri%Node1=>NodeLayerBottom(i)%NodeTri
    EdgeMid(2*i)%EdgeTri%Node2=>NodeLayerTop(i+1)%NodeTri
  END DO  
  EdgeMid(9)%EdgeTri%Node1=>NodeLayerBottom(5)%NodeTri
  EdgeMid(9)%EdgeTri%Node2=>NodeLayerTop(5)%NodeTri
  EdgeMid(10)%EdgeTri%Node1=>NodeLayerBottom(5)%NodeTri
  EdgeMid(10)%EdgeTri%Node2=>NodeLayerTop(1)%NodeTri
  DO i=1,4
    ALLOCATE(EdgeLayerBottom(i)%EdgeTri)
    EdgeLayerBottom(i)%EdgeTri%Node1=>NodeLayerBottom(i)%NodeTri
    EdgeLayerBottom(i)%EdgeTri%Node2=>NodeLayerBottom(i+1)%NodeTri
    CALL InsertEdge(EdgeLayerBottom(i)%EdgeTri,TriangularGrid)
  END DO  
  ALLOCATE(EdgeLayerBottom(5)%EdgeTri)
  EdgeLayerBottom(5)%EdgeTri%Node1=>NodeLayerBottom(5)%NodeTri
  EdgeLayerBottom(5)%EdgeTri%Node2=>NodeLayerBottom(1)%NodeTri
  CALL InsertEdge(EdgeLayerBottom(5)%EdgeTri,TriangularGrid)
  DO i=1,5
    ALLOCATE(EdgeBottom(i)%EdgeTri)
    EdgeBottom(i)%EdgeTri%Node1=>NodeBottom
    EdgeBottom(i)%EdgeTri%Node2=>NodeLayerBottom(i)%NodeTri
    CALL InsertEdge(EdgeBottom(i)%EdgeTri,TriangularGrid)
  END DO  
  DO i=1,4
    ALLOCATE(CellTop(i)%CellTri)
    CellTop(i)%CellTri%OrientE1=1
    CellTop(i)%CellTri%Edge1=>EdgeLayerTop(i)%EdgeTri
    CellTop(i)%CellTri%OrientE2=1
    CellTop(i)%CellTri%Edge2=>EdgeTop(i+1)%EdgeTri
    CellTop(i)%CellTri%OrientE3=-1
    CellTop(i)%CellTri%Edge3=>EdgeTop(i)%EdgeTri
    CALL InsertCell(CellTop(i)%CellTri,TriangularGrid)
  END DO  
  ALLOCATE(CellTop(5)%CellTri)
  CellTop(5)%CellTri%OrientE1=1
  CellTop(5)%CellTri%Edge1=>EdgeLayerTop(5)%EdgeTri
  CellTop(5)%CellTri%OrientE2=1
  CellTop(5)%CellTri%Edge2=>EdgeTop(1)%EdgeTri
  CellTop(5)%CellTri%OrientE3=-1
  CellTop(5)%CellTri%Edge3=>EdgeTop(5)%EdgeTri
  CALL InsertCell(CellTop(5)%CellTri,TriangularGrid)
  DO i=1,10
    ALLOCATE(CellMid(i)%CellTri)
    CALL InsertCell(CellMid(i)%CellTri,TriangularGrid)
  END DO  
  DO i=1,4
    CellMid(2*i-1)%CellTri%OrientE1=-1
    CellMid(2*i-1)%CellTri%Edge1=>EdgeMid(2*i-1)%EdgeTri
    CellMid(2*i-1)%CellTri%OrientE2=1
    CellMid(2*i-1)%CellTri%Edge2=>EdgeMid(2*i)%EdgeTri
    CellMid(2*i-1)%CellTri%OrientE3=-1
    CellMid(2*i-1)%CellTri%Edge3=>EdgeLayerTop(i)%EdgeTri
    CellMid(2*i)%CellTri%OrientE1=1
    CellMid(2*i)%CellTri%Edge1=>EdgeMid(2*i+1)%EdgeTri
    CellMid(2*i)%CellTri%OrientE2=-1
    CellMid(2*i)%CellTri%Edge2=>EdgeMid(2*i)%EdgeTri
    CellMid(2*i)%CellTri%OrientE3=1
    CellMid(2*i)%CellTri%Edge3=>EdgeLayerBottom(i)%EdgeTri
  END DO  
  CellMid(2*5-1)%CellTri%OrientE1=-1
  CellMid(2*5-1)%CellTri%Edge1=>EdgeMid(2*5-1)%EdgeTri
  CellMid(2*5-1)%CellTri%OrientE2=1
  CellMid(2*5-1)%CellTri%Edge2=>EdgeMid(2*5)%EdgeTri
  CellMid(2*5-1)%CellTri%OrientE3=-1
  CellMid(2*5-1)%CellTri%Edge3=>EdgeLayerTop(5)%EdgeTri
  CellMid(2*5)%CellTri%OrientE1=1
  CellMid(2*5)%CellTri%Edge1=>EdgeMid(1)%EdgeTri
  CellMid(2*5)%CellTri%OrientE2=-1
  CellMid(2*5)%CellTri%Edge2=>EdgeMid(2*5)%EdgeTri
  CellMid(2*5)%CellTri%OrientE3=1
  CellMid(2*5)%CellTri%Edge3=>EdgeLayerBottom(5)%EdgeTri
  DO i=1,4
    ALLOCATE(CellBottom(i)%CellTri)
    CellBottom(i)%CellTri%OrientE1=1
    CellBottom(i)%CellTri%Edge1=>EdgeLayerBottom(i)%EdgeTri
    CellBottom(i)%CellTri%OrientE2=-1
    CellBottom(i)%CellTri%Edge2=>EdgeBottom(i+1)%EdgeTri
    CellBottom(i)%CellTri%OrientE3=1
    CellBottom(i)%CellTri%Edge3=>EdgeBottom(i)%EdgeTri
    CALL InsertCell(CellBottom(i)%CellTri,TriangularGrid)
  END DO  
  ALLOCATE(CellBottom(5)%CellTri)
  CellBottom(5)%CellTri%OrientE1=1
  CellBottom(5)%CellTri%Edge1=>EdgeLayerBottom(5)%EdgeTri
  CellBottom(5)%CellTri%OrientE2=-1
  CellBottom(5)%CellTri%Edge2=>EdgeBottom(1)%EdgeTri
  CellBottom(5)%CellTri%OrientE3=1
  CellBottom(5)%CellTri%Edge3=>EdgeBottom(5)%EdgeTri
  CALL InsertCell(CellBottom(5)%CellTri,TriangularGrid)

END SUBROUTINE IcosahedronGrid

SUBROUTINE RefineEdge(EdgeTri)
  TYPE(EdgeTri_T), POINTER :: EdgeTri

  TYPE(EdgeTri_T), POINTER :: EdgeNew
  TYPE(NodeTri_T), POINTER :: NodeE

  NodeE=>.INS.EdgeTri%Node1
! NodeE%P=0.5d0*(EdgeTri%Node1%P+EdgeTri%Node2%P)
  NodeE%P=MidPointEdge(EdgeTri%Node1%P,EdgeTri%Node2%P)
  EdgeNew=>.INS.EdgeTri
  EdgeNew%Node1=>NodeE
  EdgeNew%Node2=>EdgeTri%Node2
  EdgeTri%Node2=>NodeE

END SUBROUTINE RefineEdge

FUNCTION MidPointEdge(P1,P2)
  TYPE(Point_T) :: MidPointEdge
  TYPE(Point_T) :: P1,P2
  REAL(8) :: Rad
  Rad=Norm(P1)
  MidPointEdge=0.5d0*(P1+P2)
  MidPointEdge=Rad/Norm(MidPointEdge)*MidPointEdge
END FUNCTION MidPointEdge

FUNCTION LastEdge(TriangularGrid)
  TYPE(EdgeTri_T), POINTER :: LastEdge
  TYPE(TriangularGrid_T) :: TriangularGrid

  LastEdge=>TriangularGrid%LastEdge
  DO
    IF (ASSOCIATED(LastEdge%Next)) THEN
      LastEdge=>LastEdge%Next
    ELSE
      EXIT
    END IF  
  END DO
END FUNCTION LastEdge

FUNCTION CountNode(TriangularGrid)
  INTEGER :: CountNode
  TYPE(TriangularGrid_T) :: TriangularGrid

  TYPE(NodeTri_T), POINTER :: StartNode

  StartNode=>TriangularGrid%StartNode
  CountNode=1
  DO
    IF (ASSOCIATED(StartNode%Next)) THEN
      StartNode=>StartNode%Next
      CountNode=CountNode+1
    ELSE
      EXIT
    END IF  
  END DO
END FUNCTION CountNode

FUNCTION CountEdge(TriangularGrid)
  INTEGER :: CountEdge
  TYPE(TriangularGrid_T) :: TriangularGrid

  TYPE(EdgeTri_T), POINTER :: StartEdge

  StartEdge=>TriangularGrid%StartEdge
  CountEdge=1
  DO
    IF (ASSOCIATED(StartEdge%Next)) THEN
      StartEdge=>StartEdge%Next
      CountEdge=CountEdge+1
    ELSE
      EXIT
    END IF  
  END DO
END FUNCTION CountEdge

FUNCTION CountCell(TriangularGrid)
  INTEGER :: CountCell
  TYPE(TriangularGrid_T) :: TriangularGrid

  TYPE(CellTri_T), POINTER :: StartCell

  StartCell=>TriangularGrid%StartCell
  CountCell=1
  DO
    IF (ASSOCIATED(StartCell%Next)) THEN
      StartCell=>StartCell%Next
      CountCell=CountCell+1
    ELSE
      EXIT
    END IF  
  END DO
END FUNCTION CountCell

SUBROUTINE WriteNode(NodeTri)
  TYPE(NodeTri_T), POINTER :: NodeTri
  WRITE(*,*) 'Node',NodeTri%Number
  WRITE(*,*) NodeTri%P
END SUBROUTINE WriteNode

SUBROUTINE WriteNodeTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(NodeTri_T), POINTER :: CurrentNode

  CurrentNode=>TriangularGrid%StartNode
  DO  
    CALL WriteNode(CurrentNode)
    IF (ASSOCIATED(CurrentNode%Next)) THEN
      CurrentNode=>CurrentNode%Next
    ELSE
      EXIT
    END IF  
  END DO
END SUBROUTINE WriteNodeTriangularGrid

SUBROUTINE WriteEdge(EdgeTri)
  TYPE(EdgeTri_T), POINTER :: EdgeTri
  WRITE(*,*) 'Edge',EdgeTri%Number
  WRITE(*,*) EdgeTri%Node1%Number
  WRITE(*,*) EdgeTri%Node2%Number
END SUBROUTINE WriteEdge

SUBROUTINE WriteEdgeTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(EdgeTri_T), POINTER :: CurrentEdge

  CurrentEdge=>TriangularGrid%StartEdge
  DO  
    CALL WriteEdge(CurrentEdge)
    IF (ASSOCIATED(CurrentEdge%Next)) THEN
      CurrentEdge=>CurrentEdge%Next
    ELSE
      EXIT
    END IF  
  END DO
END SUBROUTINE WriteEdgeTriangularGrid

SUBROUTINE WriteCell(CellTri)
  TYPE(CellTri_T), POINTER :: CellTri
  WRITE(*,*) 'Cell',CellTri%Number
  WRITE(*,*) CellTri%Edge1%Number,CellTri%Edge2%Number,CellTri%Edge3%Number
  IF (CellTri%OrientE1==1) THEN
    WRITE(*,*) CellTri%Edge1%Node1%Number
  ELSE  
    WRITE(*,*) CellTri%Edge1%Node2%Number
  END IF  
  IF (CellTri%OrientE2==1) THEN
    WRITE(*,*) CellTri%Edge2%Node1%Number
  ELSE  
    WRITE(*,*) CellTri%Edge2%Node2%Number
  END IF  
  IF (CellTri%OrientE3==1) THEN
    WRITE(*,*) CellTri%Edge3%Node1%Number
  ELSE  
    WRITE(*,*) CellTri%Edge3%Node2%Number
  END IF  
! WRITE(*,*) 'Edge1'
! WRITE(*,*) CellTri%OrientE1
! WRITE(*,*) CellTri%Edge1%Node1%P
! WRITE(*,*) CellTri%Edge1%Node2%P
! WRITE(*,*) 'Edge2'
! WRITE(*,*) CellTri%OrientE2
! WRITE(*,*) CellTri%Edge2%Node1%P
! WRITE(*,*) CellTri%Edge2%Node2%P
! WRITE(*,*) 'Edge3'
! WRITE(*,*) CellTri%OrientE3
! WRITE(*,*) CellTri%Edge3%Node1%P
! WRITE(*,*) CellTri%Edge3%Node2%P
END SUBROUTINE WriteCell

SUBROUTINE WriteCellTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(CellTri_T), POINTER :: CurrentCell

  CurrentCell=>TriangularGrid%StartCell
  DO  
    CALL WriteCell(CurrentCell)
    IF (ASSOCIATED(CurrentCell%Next)) THEN
      CurrentCell=>CurrentCell%Next
    ELSE
      EXIT
    END IF  
  END DO
END SUBROUTINE WriteCellTriangularGrid

SUBROUTINE RefineEdgeTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(EdgeTri_T), POINTER :: CurrentEdge

  CurrentEdge=>TriangularGrid%StartEdge
  DO  
    CALL RefineEdge(CurrentEdge)
    CurrentEdge=>CurrentEdge%Next
    IF (ASSOCIATED(CurrentEdge%Next)) THEN
      CurrentEdge=>CurrentEdge%Next
    ELSE
      EXIT
    END IF  
  END DO
END SUBROUTINE RefineEdgeTriangularGrid

SUBROUTINE DeallocateTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(NodeTri_T), POINTER :: StartNode,TmpNode
  TYPE(EdgeTri_T), POINTER :: StartEdge,TmpEdge
  TYPE(CellTri_T), POINTER :: StartCell,TmpCell

  StartNode=>TriangularGrid%StartNode
  DO  
    IF (ASSOCIATED(StartNode)) THEN
      DEALLOCATE(StartNode%Edge)
      TmpNode=>StartNode
      StartNode=>StartNode%Next
      DEALLOCATE(TmpNode)
    ELSE
      EXIT
    END IF  
  END DO
  StartEdge=>TriangularGrid%StartEdge
  DO  
    IF (ASSOCIATED(StartEdge)) THEN
      DEALLOCATE(StartEdge%Cell)
      TmpEdge=>StartEdge
      StartEdge=>StartEdge%Next
      DEALLOCATE(TmpEdge)
    ELSE
      EXIT
    END IF  
  END DO
  StartCell=>TriangularGrid%StartCell
  DO  
    IF (ASSOCIATED(StartCell)) THEN
      TmpCell=>StartCell
      StartCell=>StartCell%Next
      DEALLOCATE(TmpCell)
    ELSE
      EXIT
    END IF  
  END DO
END SUBROUTINE DeallocateTriangularGrid

SUBROUTINE NumberingTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(NodeTri_T), POINTER :: StartNode
  TYPE(EdgeTri_T), POINTER :: StartEdge
  TYPE(CellTri_T), POINTER :: StartCell

  INTEGER :: NodeNumber
  INTEGER :: EdgeNumber
  INTEGER :: CellNumber
  INTEGER, ALLOCATABLE :: EdgeNode(:)
  INTEGER, ALLOCATABLE :: CellEdge(:)

  StartNode=>TriangularGrid%StartNode
  NodeNumber=0
  DO  
    IF (ASSOCIATED(StartNode)) THEN
      NodeNumber=NodeNumber+1
      StartNode%Number=NodeNumber
      StartNode=>StartNode%Next
    ELSE
      EXIT
    END IF  
  END DO
  ALLOCATE(EdgeNode(NodeNumber))
  EdgeNode=0

  StartEdge=>TriangularGrid%StartEdge
  EdgeNumber=0
  DO  
    IF (ASSOCIATED(StartEdge)) THEN
      EdgeNumber=EdgeNumber+1
      StartEdge%Number=EdgeNumber
      EdgeNode(StartEdge%Node1%Number)=EdgeNode(StartEdge%Node1%Number)+1
      EdgeNode(StartEdge%Node2%Number)=EdgeNode(StartEdge%Node2%Number)+1
      StartEdge=>StartEdge%Next
    ELSE
      EXIT
    END IF  
  END DO
  ALLOCATE(CellEdge(EdgeNumber))
  CellEdge=0

  StartCell=>TriangularGrid%StartCell
  CellNumber=0
  DO  
    IF (ASSOCIATED(StartCell)) THEN
      CellNumber=CellNumber+1
      StartCell%Number=CellNumber
      CellEdge(StartCell%Edge1%Number)=CellEdge(StartCell%Edge1%Number)+1
      CellEdge(StartCell%Edge2%Number)=CellEdge(StartCell%Edge2%Number)+1
      CellEdge(StartCell%Edge3%Number)=CellEdge(StartCell%Edge3%Number)+1
      StartCell=>StartCell%Next
    ELSE
      EXIT
    END IF  
  END DO
  
  StartNode=>TriangularGrid%StartNode
  DO  
    IF (ASSOCIATED(StartNode)) THEN
      ALLOCATE(StartNode%Edge(EdgeNode(StartNode%Number)))
      StartNode=>StartNode%Next
    ELSE
      EXIT
    END IF  
  END DO
  
  StartEdge=>TriangularGrid%StartEdge
  EdgeNode=0
  DO  
    IF (ASSOCIATED(StartEdge)) THEN
      ALLOCATE(StartEdge%Cell(CellEdge(StartEdge%Number)))
      EdgeNode(StartEdge%Node1%Number)=EdgeNode(StartEdge%Node1%Number)+1
      StartEdge%Node1%Edge(EdgeNode(StartEdge%Node1%Number))=StartEdge%Number
      EdgeNode(StartEdge%Node2%Number)=EdgeNode(StartEdge%Node2%Number)+1
      StartEdge%Node2%Edge(EdgeNode(StartEdge%Node2%Number))=StartEdge%Number
      StartEdge=>StartEdge%Next
    ELSE
      EXIT
    END IF  
  END DO

  Celledge=0
  StartCell=>TriangularGrid%StartCell
  DO  
    IF (ASSOCIATED(StartCell)) THEN
      CellEdge(StartCell%Edge1%Number)=CellEdge(StartCell%Edge1%Number)+1
      StartCell%Edge1%Cell(CellEdge(StartCell%Edge1%Number))=StartCell%Number
      CellEdge(StartCell%Edge2%Number)=CellEdge(StartCell%Edge2%Number)+1
      StartCell%Edge2%Cell(CellEdge(StartCell%Edge2%Number))=StartCell%Number
      CellEdge(StartCell%Edge3%Number)=CellEdge(StartCell%Edge3%Number)+1
      StartCell%Edge3%Cell(CellEdge(StartCell%Edge3%Number))=StartCell%Number
      StartCell=>StartCell%Next
    ELSE
      EXIT
    END IF  
  END DO
  DEALLOCATE(EdgeNode)
  DEALLOCATE(CellEdge)

END SUBROUTINE NumberingTriangularGrid

FUNCTION MidPoint(CellTri)
  TYPE(Point_T) :: MidPoint
  TYPE(CellTri_T), POINTER :: CellTri

  MidPoint=0.0d0
  IF (CellTri%OrientE1==1) THEN
    MidPoint=CellTri%Edge1%Node1%P
  ELSE  
    MidPoint=CellTri%Edge1%Node2%P
  END IF  
  IF (CellTri%OrientE2==1) THEN
    MidPoint=MidPoint+CellTri%Edge2%Node1%P
  ELSE  
    MidPoint=MidPoint+CellTri%Edge2%Node2%P
  END IF  
  IF (CellTri%OrientE3==1) THEN
    MidPoint=MidPoint+CellTri%Edge3%Node1%P
  ELSE  
    MidPoint=MidPoint+CellTri%Edge3%Node2%P
  END IF  
  MidPoint=(1.0d0/3.0d0)*MidPoint
END FUNCTION MidPoint

SUBROUTINE RefineCellTriangularGrid(TriangularGrid)
  TYPE(TriangularGrid_T) :: TriangularGrid
  TYPE(CellTri_T), POINTER :: CurrentCell

  CurrentCell=>TriangularGrid%StartCell
  DO  
    CALL RefineCell(CurrentCell)
    CurrentCell=>CurrentCell%Next
    CurrentCell=>CurrentCell%Next
    CurrentCell=>CurrentCell%Next
    IF (ASSOCIATED(CurrentCell%Next)) THEN
      CurrentCell=>CurrentCell%Next
    ELSE
      EXIT
    END IF  
  END DO
END SUBROUTINE RefineCellTriangularGrid

SUBROUTINE InsertCell(CellTri,TriangularGrid)
  TYPE(CellTri_T), POINTER :: CellTri
  TYPE(TriangularGrid_T) :: TriangularGrid

  IF (.NOT.ASSOCIATED(TriangularGrid%StartCell)) THEN
    TriangularGrid%StartCell=>CellTri
    TriangularGrid%CurrentCell=>CellTri
  ELSE  
    TriangularGrid%CurrentCell%Next=>CellTri
    TriangularGrid%CurrentCell=>CellTri
  END IF  
END SUBROUTINE InsertCell

SUBROUTINE InsertEdge(EdgeTri,TriangularGrid)
  TYPE(EdgeTri_T), POINTER :: EdgeTri
  TYPE(TriangularGrid_T) :: TriangularGrid

  IF (.NOT.ASSOCIATED(TriangularGrid%StartEdge)) THEN
    TriangularGrid%StartEdge=>EdgeTri
    TriangularGrid%CurrentEdge=>EdgeTri
  ELSE  
    TriangularGrid%CurrentEdge%Next=>EdgeTri
    TriangularGrid%CurrentEdge=>EdgeTri
  END IF  
END SUBROUTINE InsertEdge

SUBROUTINE InsertNode(NodeTri,TriangularGrid)
  TYPE(NodeTri_T), POINTER :: NodeTri
  TYPE(TriangularGrid_T) :: TriangularGrid

  IF (.NOT.ASSOCIATED(TriangularGrid%StartNode)) THEN
    TriangularGrid%StartNode=>NodeTri
    TriangularGrid%CurrentNode=>NodeTri
  ELSE  
    TriangularGrid%CurrentNode%Next=>NodeTri
    TriangularGrid%CurrentNode=>NodeTri
  END IF  
END SUBROUTINE InsertNode

FUNCTION AppendNode(Node)
  TYPE(NodeTri_T), POINTER :: AppendNode
  TYPE(NodeTri_T), POINTER, INTENT(In) :: Node

  ALLOCATE(AppendNode)
  AppendNode%Next=>Node%Next
  Node%Next=>AppendNode
END FUNCTION AppendNode

FUNCTION AppendEdge(Edge)
  TYPE(EdgeTri_T), POINTER :: AppendEdge
  TYPE(EdgeTri_T), POINTER, INTENT(In) :: Edge

  ALLOCATE(AppendEdge)
  AppendEdge%Next=>Edge%Next
  Edge%Next=>AppendEdge
END FUNCTION AppendEdge

FUNCTION AppendCell(Cell)
  TYPE(CellTri_T), POINTER :: AppendCell
  TYPE(CellTri_T), POINTER, INTENT(In) :: Cell

  ALLOCATE(AppendCell)
  AppendCell%Next=>Cell%Next
  Cell%Next=>AppendCell
END FUNCTION AppendCell

SUBROUTINE RefineCell(CellTri)
 TYPE(CellTri_T), POINTER :: CellTri

 TYPE(EdgeTri_T), POINTER :: Edge1,Edge2,Edge3
 TYPE(EdgeTri_T), POINTER :: Edge1N,Edge2N,Edge3N
 TYPE(EdgeTri_T), POINTER :: EdgeI1,EdgeI2,EdgeI3
 TYPE(CellTri_T), POINTER :: Cell1,Cell2,Cell3

 Edge1=>CellTri%Edge1
 Edge1N=>CellTri%Edge1%Next
 Edge2=>CellTri%Edge2
 Edge2N=>CellTri%Edge2%Next
 Edge3=>CellTri%Edge3
 Edge3N=>CellTri%Edge3%Next

! Insert new cells 
 Cell1=>.INS.CellTri
 Cell2=>.INS.CellTri
 Cell3=>.INS.CellTri
  
 IF (CellTri%OrientE1==1) THEN 
   Cell2%Edge2=>Edge1 
   Cell2%OrientE2=1
   Cell3%Edge1=>Edge1N 
   Cell3%OrientE1=1
 ELSE  
   Cell2%Edge2=>Edge1N 
   Cell2%OrientE2=-1
   Cell3%Edge1=>Edge1 
   Cell3%OrientE1=-1
 END IF  
 IF (CellTri%OrientE2==1) THEN 
   Cell3%Edge2=>Edge2 
   Cell3%OrientE2=1
   Cell1%Edge1=>Edge2N 
   Cell1%OrientE1=1
 ELSE  
   Cell3%Edge2=>Edge2N 
   Cell3%OrientE2=-1
   Cell1%Edge1=>Edge2 
   Cell1%OrientE1=-1
 END IF  
 IF (CellTri%OrientE3==1) THEN 
   Cell1%Edge2=>Edge3 
   Cell1%OrientE2=1
   Cell2%Edge1=>Edge3N 
   Cell2%OrientE1=1
 ELSE  
   Cell1%Edge2=>Edge3N 
   Cell1%OrientE2=-1
   Cell2%Edge1=>Edge3 
   Cell2%OrientE1=-1
 END IF  
 EdgeI1=>.INS.Edge1N
 EdgeI2=>.INS.Edge2N
 EdgeI3=>.INS.Edge3N
 EdgeI1%Node1=>Edge2%Node2
 EdgeI1%Node2=>Edge3%Node2
 EdgeI2%Node1=>Edge3%Node2
 EdgeI2%Node2=>Edge1%Node2
 EdgeI3%Node1=>Edge1%Node2
 EdgeI3%Node2=>Edge2%Node2
 Cell1%Edge3=>EdgeI1
 Cell1%OrientE3=-1
 Cell2%Edge3=>EdgeI2
 Cell2%OrientE3=-1
 Cell3%Edge3=>EdgeI3
 Cell3%OrientE3=-1
 CellTri%Edge1=>EdgeI1
 CellTri%OrientE1=1
 CellTri%Edge2=>EdgeI2
 CellTri%OrientE2=1
 CellTri%Edge3=>EdgeI3
 CellTri%OrientE3=1
END SUBROUTINE RefineCell

END MODULE Triangular_Mod
=#
