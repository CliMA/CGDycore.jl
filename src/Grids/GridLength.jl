function GridLength(Grid,::Grids.Quad) 
  Faces = Grid.Faces
  LengthLoc = 0
  for iF = 1 : Grid.NumFaces
    LengthLoc = max(sqrt(Faces[iF].Area), LengthLoc)
  end
  Length = MPI.Allreduce(LengthLoc, MPI.MAX, MPI.COMM_WORLD)
  return Length
end
