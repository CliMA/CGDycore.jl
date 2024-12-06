function GridLength(Grid) 
  Faces = Grid.Faces
  LengthMaxLoc = 0
  LengthMinLoc = 4 * pi * Grid.Rad^2
  for iF = 1 : Grid.NumFaces
    LengthMaxLoc = max(sqrt(Faces[iF].Area), LengthMaxLoc)
    LengthMinLoc = min(sqrt(Faces[iF].Area), LengthMinLoc)
  end
  LengthMax = MPI.Allreduce(LengthMaxLoc, MPI.MAX, MPI.COMM_WORLD)
  LengthMin = MPI.Allreduce(LengthMinLoc, MPI.MIN, MPI.COMM_WORLD)
  return LengthMin,LengthMax
end
