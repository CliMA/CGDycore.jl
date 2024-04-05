function ProjectFace!(backend,FTB,p,Grid,F)
  x = zeros(3)
  for iF = 1 : Grid.NumFaces
    x[1] = Grid.Faces[iF].Mid.x  
    x[2] = Grid.Faces[iF].Mid.y  
    x[3] = Grid.Faces[iF].Mid.z  
    p[iF], = F(x,0.0)
  end
end
