function ProjectFace!(backend,FTB,p,Grid,F)
  x = zeros(3)
  for iF = 1 : Grid.NumFaces
    x[1] = Grid.Faces[iF].Mid.x  
    x[2] = Grid.Faces[iF].Mid.y  
    x[3] = Grid.Faces[iF].Mid.z  
    p[iF], = F(x,0.0)
  end
end

function ProjectFaceTr!(backend,FTB,p,Grid,F)
  x = zeros(3)
  for iF = 1 : Grid.NumFaces
    x[1] = Grid.Faces[iF].Mid.x
    x[2] = Grid.Faces[iF].Mid.y
    x[3] = Grid.Faces[iF].Mid.z
    _,_,_,_,p[iF] = F(x,0.0)
  end
end

function ProjectEdgeScalar!(backend,FTB,p,Grid,F)
  x = zeros(3)
  for iF = 1 : Grid.NumEdges
    x[1] = Grid.Edges[iF].Mid.x
    x[2] = Grid.Edges[iF].Mid.y
    x[3] = Grid.Edges[iF].Mid.z
    p[iF], = F(x,0.0)
  end
end

function ProjectEdge!(backend,FTB,u,Grid,F)
  x = zeros(3)
  n = zeros(3)
  VelSp = zeros(3)  
  for iE = 1 : Grid.NumEdges
    x[1] = Grid.Edges[iE].Mid.x
    x[2] = Grid.Edges[iE].Mid.y
    x[3] = Grid.Edges[iE].Mid.z
    _,VelSp[1],VelSp[2],VelSp[3], = F(x,0.0)
    lon,lat,r = Grids.cart2sphere(x[1],x[2],x[3])
    VelCa = VelSphere2Cart(VelSp,lon,lat)
    n[1] = Grid.Edges[iE].n.x
    n[2] = Grid.Edges[iE].n.y
    n[3] = Grid.Edges[iE].n.z
    u[iE] = n[1] * VelCa[1] + n[2] * VelCa[2] + n[3] * VelCa[3]
  end
end

function ProjectDualEdge!(backend,FTB,u,Grid,F)
  x = zeros(3)
  n = zeros(3)
  VelSp = zeros(3)
  for iE = 1 : Grid.NumEdges
    x[1] = Grid.Edges[iE].Mid.x
    x[2] = Grid.Edges[iE].Mid.y
    x[3] = Grid.Edges[iE].Mid.z
    _,VelSp[1],VelSp[2],VelSp[3], = F(x,0.0)
    lon,lat,r = Grids.cart2sphere(x[1],x[2],x[3])
    VelCa = VelSphere2Cart(VelSp,lon,lat)
    u[iE] = Grid.Edges[iE].t.x * VelCa[1] + Grid.Edges[iE].t.y * VelCa[2] + Grid.Edges[iE].t.z * VelCa[3]
  end
end

function ConvertVelocityCart!(backend,FTB,VelCa,Vel,Grid)
  for iF = 1 : Grid.NumFaces
    N = zeros(length(Grid.Faces[iF].E)+1,3)  
    Rhs = zeros(length(Grid.Faces[iF].E)+1)
    k = zeros(3)
    t = zeros(3)
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      N[i,1] = Grid.Edges[iE].n.x
      N[i,2] = Grid.Edges[iE].n.y
      N[i,3] = Grid.Edges[iE].n.z
      Rhs[i] = Vel[iE]
    end
    k[1] = Grid.Faces[iF].Mid.x
    k[2] = Grid.Faces[iF].Mid.y
    k[3] = Grid.Faces[iF].Mid.z
    k = k / norm(k)
    N[length(Grid.Faces[iF].E)+1,:] = k
    Rhs[length(Grid.Faces[iF].E)+1] = 0
    VelCa[iF,:] = N \ Rhs
  end
end

function ConvertVelocityTCart!(backend,FTB,VelCa,VelT,Grid)
  for iF = 1 : Grid.NumFaces
    T = zeros(length(Grid.Faces[iF].E)+1,3)
    Rhs = zeros(length(Grid.Faces[iF].E)+1)
    k = zeros(3)
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      T[i,1] = Grid.Edges[iE].t.x
      T[i,2] = Grid.Edges[iE].t.y
      T[i,3] = Grid.Edges[iE].t.z
      Rhs[i] = VelT[iE]
    end
    k[1] = Grid.Faces[iF].Mid.x
    k[2] = Grid.Faces[iF].Mid.y
    k[3] = Grid.Faces[iF].Mid.z
    k = k / norm(k)
    T[length(Grid.Faces[iF].E)+1,:] = k
    Rhs[length(Grid.Faces[iF].E)+1] = 0
    VelCa[iF,:] = T \ Rhs
  end
end

function ConvertVelocitySp!(backend,FTB,VelSp,Vel,Grid)
  for iF = 1 : Grid.NumFaces
    N = zeros(length(Grid.Faces[iF].E)+1,3)
    Rhs = zeros(length(Grid.Faces[iF].E)+1)
    k = zeros(3)
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      N[i,1] = Grid.Edges[iE].n.x
      N[i,2] = Grid.Edges[iE].n.y
      N[i,3] = Grid.Edges[iE].n.z
      Rhs[i] = Vel[iE]
    end
    k[1] = Grid.Faces[iF].Mid.x
    k[2] = Grid.Faces[iF].Mid.y
    k[3] = Grid.Faces[iF].Mid.z
    k = k / norm(k)
    N[length(Grid.Faces[iF].E)+1,:] = k
    Rhs[length(Grid.Faces[iF].E)+1] = 0
    VelCa = N \ Rhs
    lon,lat,_ = Grids.cart2sphere(Grid.Faces[iF].Mid.x,Grid.Faces[iF].Mid.y,Grid.Faces[iF].Mid.z)
    VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end
end

function ConvertVelocityTSp!(backend,FTB,VelSp,VelT,Grid)
  for iF = 1 : Grid.NumFaces
    T = zeros(length(Grid.Faces[iF].E)+1,3)
    Rhs = zeros(length(Grid.Faces[iF].E)+1)
    k = zeros(3)
    for i = 1 : length(Grid.Faces[iF].E)
      iE = Grid.Faces[iF].E[i]
      T[i,1] = Grid.Edges[iE].t.x
      T[i,2] = Grid.Edges[iE].t.y
      T[i,3] = Grid.Edges[iE].t.z
      Rhs[i] = VelT[iE]
    end
    k[1] = Grid.Faces[iF].Mid.x
    k[2] = Grid.Faces[iF].Mid.y
    k[3] = Grid.Faces[iF].Mid.z
    k = k / norm(k)
    T[length(Grid.Faces[iF].E)+1,:] = k
    Rhs[length(Grid.Faces[iF].E)+1] = 0
    VelCa = T \ Rhs
    lon,lat,_ = Grids.cart2sphere(Grid.Faces[iF].Mid.x,Grid.Faces[iF].Mid.y,Grid.Faces[iF].Mid.z)
    VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
    VelSp[iF,1] = VelSpLoc[1]
    VelSp[iF,2] = VelSpLoc[2]
  end
end

function VelCart2Sphere(VelCa,lon,lat)

  Rotate = zeros(3,3)
  Rotate[1,1] =           -sin(lon)
  Rotate[2,1] =-sin(lat) * cos(lon)
  Rotate[3,1] = cos(lat) * cos(lon)
  Rotate[1,2] =           cos(lon)
  Rotate[2,2] =-sin(lat) * sin(lon)
  Rotate[3,2] = cos(lat) * sin(lon)
  Rotate[1,3] =          0.0
  Rotate[2,3] = cos(lat)
  Rotate[3,3] = sin(lat)
  VelSp = zeros(3,1)
  VelSp[1] = Rotate[1,1] * VelCa[1] +
    Rotate[1,2] * VelCa[2] +
    Rotate[1,3] * VelCa[3]
  VelSp[2] = Rotate[2,1] * VelCa[1] +
    Rotate[2,2] * VelCa[2] +
    Rotate[2,3] * VelCa[3]
  VelSp[3] = Rotate[3,1]*VelCa[1] +
    Rotate[3,2] * VelCa[2] +
    Rotate[3,3] * VelCa[3]
  return VelSp
end

function VelSphere2Cart(VelSp,lon,lat)
  rot=zeros(3,3);
  rot[1,1]=-sin(lon);
  rot[2,1]=-sin(lat)*cos(lon);
  rot[3,1]= cos(lat)*cos(lon);
  rot[1,2]=cos(lon);
  rot[2,2]=-sin(lat)*sin(lon);
  rot[3,2]=cos(lat)*sin(lon);
  rot[1,3]=0.0;
  rot[2,3]=cos(lat);
  rot[3,3]=sin(lat);
  VelCa = rot'*VelSp
end
