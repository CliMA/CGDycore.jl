function ProjectFace(backend,FTB,Grid,F)
  p=zeros(Grid.NumFaces)
  for iF = 1 : Grid.NumFaces
    x = Grid.Faces[iF].Mid.x  
    y = Grid.Faces[iF].Mid.y  
    z = Grid.Faces[iF].Mid.z  
    p[iF] = F(x,y,z)
  end
  return p
end
