function ProjectFace!(backend,FTB,p,Grid,F)
  x = zeros(3)
  for iF = 1 : Grid.NumFaces
    x[1] = Grid.Faces[iF].Mid.x  
    x[2] = Grid.Faces[iF].Mid.y  
    x[3] = Grid.Faces[iF].Mid.z  
    p[iF], = F(x,0.0)
  end
end

function ProjectCell3D!(backend,FTB,p,Metric,Grid,F)
  x = zeros(3)
  for iF = 1 : Grid.NumFaces
    for iz = 1 : Grid.nz  
      x[1] = Metric.PrimalMidPoints[1,iz,iF]
      x[2] = Metric.PrimalMidPoints[2,iz,iF]
      x[3] = Metric.PrimalMidPoints[3,iz,iF]
      p[iz,iF], = F(x,0.0)
    end
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

function ProjectFace3D!(backend,FTB,u,Metric,Grid,F)
  x = zeros(3)
  n = zeros(3)
  VelSp = zeros(3)
  uS,uT = u
  for iE = 1 : Grid.NumEdges
    for iz = 1 : Grid.nz   
      x[1] = Metric.PrimalSideMidPoints[1,iz,iE]  
      x[2] = Metric.PrimalSideMidPoints[2,iz,iE]  
      x[3] = Metric.PrimalSideMidPoints[3,iz,iE]  
      _,VelSp[1],VelSp[2],VelSp[3], = F(x,0.0)
      lon,lat,r = Grids.cart2sphere(x[1],x[2],x[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      n[1] = Metric.PrimalSideNormals[1,iz,iE]
      n[2] = Metric.PrimalSideNormals[2,iz,iE]
      n[3] = Metric.PrimalSideNormals[3,iz,iE]
      uS[iz,iE] = n[1] * VelCa[1] + n[2] * VelCa[2] + n[3] * VelCa[3]
    end
  end
  for iF = 1 : Grid.NumFaces
    @. @views uT[1,:] = 0.0  
    @. @views uT[Grid.nz+1,:] = 0.0  
    for iz = 2 : Grid.nz   
      x[1] = Metric.PrimalTopMidPoints[1,iz,iF]  
      x[2] = Metric.PrimalTopMidPoints[2,iz,iF]  
      x[2] = Metric.PrimalTopMidPoints[3,iz,iF]  
      _,VelSp[1],VelSp[2],VelSp[3], = F(x,0.0)
      lon,lat,r = Grids.cart2sphere(x[1],x[2],x[3])
      VelCa = VelSphere2Cart(VelSp,lon,lat)
      n[1] = Metric.PrimalTopNormals[1,iz,iF]  
      n[2] = Metric.PrimalTopNormals[2,iz,iF]  
      n[3] = Metric.PrimalTopNormals[3,iz,iF]  
      uT[iz,iF] = n[1] * VelCa[1] + n[2] * VelCa[2] + n[3] * VelCa[3]
    end
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

function ConvertVelocity3DSp!(backend,FTB,VelSp,Vel,Metric,Grid)
  N = zeros(6,3)
  Rhs = zeros(6)
  VelS,VelT = Vel
  iSp = 0
  for iF = 1 : Grid.NumFaces
    for iz = 1 : Grid.nz  
      iE1 = Grid.Faces[iF].E[1]   
      @. @views N[1,:] = Metric.PrimalSideNormals[:,iz,iE1]
      Rhs[1] = VelS[iz,iE1]
      iE2 = Grid.Faces[iF].E[2]   
      @. @views N[2,:] = Metric.PrimalSideNormals[:,iz,iE2]
      Rhs[2] = VelS[iz,iE2]
      iE3 = Grid.Faces[iF].E[3]   
      @. @views N[3,:] = Metric.PrimalSideNormals[:,iz,iE3]
      Rhs[3] = VelS[iz,iE3]
      iE4 = Grid.Faces[iF].E[4]   
      @. @views N[4,:] = Metric.PrimalSideNormals[:,iz,iE4]
      Rhs[4] = VelS[iz,iE4]
      @. @views N[5,:] = Metric.PrimalTopNormals[:,iz,iF]
      Rhs[5] = VelT[iz,iF]
      @. @views N[6,:] = Metric.PrimalTopNormals[:,iz+1,iF]
      Rhs[6] = VelT[iz+1,iF]
      VelCa = N \ Rhs
      x1 = Metric.PrimalMidPoints[1,iz,iF]
      x2 = Metric.PrimalMidPoints[2,iz,iF]
      x3 = Metric.PrimalMidPoints[3,iz,iF]
      lon,lat,_ = Grids.cart2sphere(x1,x2,x3)
      VelSpLoc = VelCart2Sphere(VelCa,lon,lat)
      iSp += 1
      VelSp[iSp,1] = VelSpLoc[1]
      VelSp[iSp,2] = VelSpLoc[2]
      VelSp[iSp,3] = VelSpLoc[3]
    end
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
