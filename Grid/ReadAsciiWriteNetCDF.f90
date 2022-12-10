PROGRAM ReadAsciiWriteNetCDF

  USE netcdf

  IMPLICIT NONE

  CHARACTER(80) :: InputFilename
  INTEGER :: i,iVert,iEdge,iCell
  INTEGER :: NumVert,NumVert_id,Vertices_id
  REAL(8), ALLOCATABLE :: Vertices(:,:)
  INTEGER :: NumEdges,NumEdges_id,edges_id
  REAL(8), ALLOCATABLE :: Edges(:,:)
  INTEGER :: NumCells,NumCells_id,Cells_id
  REAL(8), ALLOCATABLE :: Cells(:,:)
  INTEGER :: info
  INTEGER :: ncid
  INTEGER :: dim,dim_id
  INTEGER :: nv,nv_id
  INTEGER :: ne,ne_id

  CALL GET_COMMAND_ARGUMENT(1,InputFileName)

  dim=2
  nv=2
  ne=4
  OPEN(UNIT=10,File=TRIM(InputFileName)//'.txt',STATUS='UNKNOWN')
  READ(10,*) NumVert
  ALLOCATE(Vertices(dim,NumVert))
  DO i=1,NumVert
    READ(10,*) iVert,Vertices(:,i) 
  END DO  
  READ(10,*) NumEdges
  ALLOCATE(Edges(nv,NumEdges))
  DO i=1,NumEdges
    READ(10,*) iEdge,Edges(:,i) 
  END DO  
  READ(10,*) NumCells
  ALLOCATE(Cells(ne,NumCells))
  DO i=1,NumCells
    READ(10,*) iCell,Cells(:,i) 
  END DO  
  CLOSE(UNIT=10)

  info=nf90_create(TRIM(InputFileName)//'.nc',NF90_CLOBBER, ncid) 
  info=nf90_def_dim(ncid, "dim", dim, dim_id) 
  info=nf90_def_dim(ncid, "nv", nv, nv_id) 
  info=nf90_def_dim(ncid, "ne", ne, ne_id) 
  info=nf90_def_dim(ncid, "NumVert", NumVert, NumVert_id) 
  info=nf90_def_dim(ncid, "NumEdges", NumEdges, NumEdges_id) 
  info=nf90_def_dim(ncid, "NumCells", NumCells, NumCells_id) 
  info=nf90_def_var(ncid, "Vertices", NF90_DOUBLE, (/dim_id,NumVert_id/), Vertices_id)
  info=nf90_def_var(ncid, "Edges", NF90_INT, (/nv_id,NumEdges_id/), Edges_id)
  info=nf90_def_var(ncid, "Cells", NF90_INT, (/ne_id,NumCells_id/), Cells_id)
  info=nf90_enddef(ncid)

  info=nf90_put_var(ncid, Vertices_id, Vertices)
  info=nf90_put_var(ncid, Edges_id, Edges)
  info=nf90_put_var(ncid, Cells_id, Cells)
  info=nf90_close(ncid)

END PROGRAM ReadAsciiWriteNetCDF
