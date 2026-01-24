#TriSphere   
@inline function Jacobi!(J,detJ,pinvJ,X,::Grids.TriPlanar,ksi1,ksi2,F,Grid)
  Rad =   Grid.Rad

  XT1 = 0.5 * (Grid.Nodes[F.N[1]].P.x*(-ksi1 - ksi2)+
        Grid.Nodes[F.N[2]].P.x * (1 + ksi1) +
        Grid.Nodes[F.N[3]].P.x * (1 + ksi2))

  XT2 = 0.5 * (Grid.Nodes[F.N[1]].P.y*(-ksi1 - ksi2)+
        Grid.Nodes[F.N[2]].P.y * (1 + ksi1) +
        Grid.Nodes[F.N[3]].P.y * (1 + ksi2))
       
  XT3 = 0.5 * (Grid.Nodes[F.N[1]].P.z * (-ksi1 - ksi2)+
        Grid.Nodes[F.N[2]].P.z * (1 + ksi1) +
        Grid.Nodes[F.N[3]].P.z * (1 + ksi2))

  XLoc = SVector{3}(XT1,XT2,XT3)

  JP = @SArray[Grid.Nodes[F.N[1]].P.x Grid.Nodes[F.N[2]].P.x Grid.Nodes[F.N[3]].P.x;
               Grid.Nodes[F.N[1]].P.y Grid.Nodes[F.N[2]].P.y Grid.Nodes[F.N[3]].P.y;
               Grid.Nodes[F.N[1]].P.z Grid.Nodes[F.N[2]].P.z Grid.Nodes[F.N[3]].P.z]

  J3 = @SArray([-1/2 -1/2;
                 1/2   0;
                  0   1/2])


  J .= JP * J3      
  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ  .= pinvJac(J)
  X .= XLoc 
end

@inline function Jacobi!(J,detJ,pinvJ,X,::Grids.Tri,ksi1,ksi2,F,Grid)
  Rad =   Grid.Rad

  XT1 = 0.5 * (Grid.Nodes[F.N[1]].P.x*(-ksi1 - ksi2)+
        Grid.Nodes[F.N[2]].P.x * (1 + ksi1) +
        Grid.Nodes[F.N[3]].P.x * (1 + ksi2))

  XT2 = 0.5 * (Grid.Nodes[F.N[1]].P.y*(-ksi1 - ksi2)+
        Grid.Nodes[F.N[2]].P.y * (1 + ksi1) +
        Grid.Nodes[F.N[3]].P.y * (1 + ksi2))
       
  XT3 = 0.5 * (Grid.Nodes[F.N[1]].P.z * (-ksi1 - ksi2)+
        Grid.Nodes[F.N[2]].P.z * (1 + ksi1) +
        Grid.Nodes[F.N[3]].P.z * (1 + ksi2))

  XLoc = SVector{3}(XT1,XT2,XT3)

  JP = @SArray[Grid.Nodes[F.N[1]].P.x Grid.Nodes[F.N[2]].P.x Grid.Nodes[F.N[3]].P.x;
               Grid.Nodes[F.N[1]].P.y Grid.Nodes[F.N[2]].P.y Grid.Nodes[F.N[3]].P.y;
               Grid.Nodes[F.N[1]].P.z Grid.Nodes[F.N[2]].P.z Grid.Nodes[F.N[3]].P.z]

  J3 = @SArray([-1/2 -1/2;
                 1/2   0;
                  0   1/2])


  f=Rad*(XT1^2+XT2^2+XT3^2)^(-3/2)
  dX1dXT1=f*(XT2^2+XT3^2)
  dX1dXT2=-f*XT1*XT2 
  dX1dXT3=-f*XT1*XT3
  dX2dXT1=dX1dXT2 
  dX2dXT2=f*(XT1^2+XT3^2)
  dX2dXT3=-f*XT2*XT3
  dX3dXT1=dX1dXT3 
  dX3dXT2=dX2dXT3 
  dX3dXT3  =f*(XT1^2+XT2^2)

  J1 = @SArray([dX1dXT1    dX2dXT1     dX3dXT1   
                dX1dXT2     dX2dXT2     dX3dXT2
                dX1dXT3     dX2dXT3     dX3dXT3])   
  J .= J1 * JP * J3      
  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ  .= pinvJac(J)
  X .= XLoc / norm(XLoc) * Rad
end

@inline function Jacobi!(J,detJ,pinvJ,X,::Grids.Tri,ksi1,ksi2,P1,P2,P3,Rad)

  XT1 = 0.5 * (P1.x*(-ksi1 - ksi2)+
        P2.x * (1 + ksi1) +
        P3.x * (1 + ksi2))

  XT2 = 0.5 * (P1.y*(-ksi1 - ksi2)+
        P2.y * (1 + ksi1) +
        P3.y * (1 + ksi2))
       
  XT3 = 0.5 * (P1.z * (-ksi1 - ksi2)+
        P2.z * (1 + ksi1) +
        P3.z * (1 + ksi2))

  XLoc = SVector{3}(XT1,XT2,XT3)

  JP = @SArray[P1.x P2.x P3.x;
               P1.y P2.y P3.y;
               P1.z P2.z P3.z]

  J3 = @SArray([-1/2 -1/2;
                 1/2   0;
                  0   1/2])


  f=Rad*(XT1^2+XT2^2+XT3^2)^(-3/2)
  dX1dXT1=f*(XT2^2+XT3^2)
  dX1dXT2=-f*XT1*XT2 
  dX1dXT3=-f*XT1*XT3
  dX2dXT1=dX1dXT2 
  dX2dXT2=f*(XT1^2+XT3^2)
  dX2dXT3=-f*XT2*XT3
  dX3dXT1=dX1dXT3 
  dX3dXT2=dX2dXT3 
  dX3dXT3  =f*(XT1^2+XT2^2)

  J1 = @SArray([dX1dXT1    dX2dXT1     dX3dXT1   
                dX1dXT2     dX2dXT2     dX3dXT2
                dX1dXT3     dX2dXT3     dX3dXT3])   
  J .= J1 * JP * J3      
  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ  .= pinvJac(J)
  X .= XLoc / norm(XLoc) * Rad
end


@inline function Jacobi!(J,detJ,pinvJ,X,::Grids.Quad,ksi1,ksi2,F,Grid)
  Rad =   Grid.Rad
  
  XT1 = 0.25*(Grid.Nodes[F.N[1]].P.x*(1-ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[2]].P.x*(1+ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[3]].P.x*(1+ksi1)*(1+ksi2)+
              Grid.Nodes[F.N[4]].P.x*(1-ksi1)*(1+ksi2))
    
  XT2 = 0.25*(Grid.Nodes[F.N[1]].P.y*(1-ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[2]].P.y*(1+ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[3]].P.y*(1+ksi1)*(1+ksi2)+
              Grid.Nodes[F.N[4]].P.y*(1-ksi1)*(1+ksi2))
           
  XT3 = 0.25*(Grid.Nodes[F.N[1]].P.z*(1-ksi1)*(1-ksi2)+
             Grid.Nodes[F.N[2]].P.z*(1+ksi1)*(1-ksi2)+
             Grid.Nodes[F.N[3]].P.z*(1+ksi1)*(1+ksi2)+
             Grid.Nodes[F.N[4]].P.z*(1-ksi1)*(1+ksi2))
    
  XLoc = SVector{3}(XT1,XT2,XT3)
  JP = @SArray[Grid.Nodes[F.N[1]].P.x Grid.Nodes[F.N[2]].P.x Grid.Nodes[F.N[3]].P.x Grid.Nodes[F.N[4]].P.x;
               Grid.Nodes[F.N[1]].P.y Grid.Nodes[F.N[2]].P.y Grid.Nodes[F.N[3]].P.y Grid.Nodes[F.N[4]].P.y;
               Grid.Nodes[F.N[1]].P.z Grid.Nodes[F.N[2]].P.z Grid.Nodes[F.N[3]].P.z Grid.Nodes[F.N[4]].P.z]

  J3 = @SArray([-0.25 + 0.25*ksi2  -0.25 + 0.25*ksi1
                 0.25 - 0.25*ksi2  -0.25 - 0.25*ksi1
                 0.25 + 0.25*ksi2   0.25 + 0.25*ksi1
                -0.25 - 0.25*ksi2   0.25 - 0.25*ksi1])


  f = Rad * (XT1^2 + XT2^2 + XT3^2)^(-3/2)
  dX1dXT1 = f * (XT2^2 + XT3^2)
  dX1dXT2 = -f * XT1 * XT2 
  dX1dXT3 = -f * XT1 * XT3
  dX2dXT1 = dX1dXT2 
  dX2dXT2 = f *(XT1^2 + XT3^2)
  dX2dXT3 = -f * XT2 * XT3
  dX3dXT1 = dX1dXT3 
  dX3dXT2 = dX2dXT3 
  dX3dXT3 = f * (XT1^2 + XT2^2)

  J1  =   @SArray([dX1dXT1    dX1dXT2     dX1dXT3   
                dX2dXT1     dX2dXT2     dX2dXT3
                dX3dXT1     dX3dXT2     dX3dXT3])   
  J   .=   J1*JP*J3      

  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ  .= pinvJac(J)
  X .= XLoc / norm(XLoc) * Rad
end

@inline function Jacobi!(J,X,::Grids.Quad,ksi1,ksi2,F,Grid)
  Rad =   Grid.Rad
  
  XT1 = 0.25*(Grid.Nodes[F.N[1]].P.x*(1-ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[2]].P.x*(1+ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[3]].P.x*(1+ksi1)*(1+ksi2)+
              Grid.Nodes[F.N[4]].P.x*(1-ksi1)*(1+ksi2))
    
  XT2 = 0.25*(Grid.Nodes[F.N[1]].P.y*(1-ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[2]].P.y*(1+ksi1)*(1-ksi2)+
              Grid.Nodes[F.N[3]].P.y*(1+ksi1)*(1+ksi2)+
              Grid.Nodes[F.N[4]].P.y*(1-ksi1)*(1+ksi2))
           
  XT3 = 0.25*(Grid.Nodes[F.N[1]].P.z*(1-ksi1)*(1-ksi2)+
             Grid.Nodes[F.N[2]].P.z*(1+ksi1)*(1-ksi2)+
             Grid.Nodes[F.N[3]].P.z*(1+ksi1)*(1+ksi2)+
             Grid.Nodes[F.N[4]].P.z*(1-ksi1)*(1+ksi2))
    
  XLoc = SVector{3}(XT1,XT2,XT3)
  JP = @SArray[Grid.Nodes[F.N[1]].P.x Grid.Nodes[F.N[2]].P.x Grid.Nodes[F.N[3]].P.x Grid.Nodes[F.N[4]].P.x;
               Grid.Nodes[F.N[1]].P.y Grid.Nodes[F.N[2]].P.y Grid.Nodes[F.N[3]].P.y Grid.Nodes[F.N[4]].P.y;
               Grid.Nodes[F.N[1]].P.z Grid.Nodes[F.N[2]].P.z Grid.Nodes[F.N[3]].P.z Grid.Nodes[F.N[4]].P.z]

  J3 = @SArray([-1 + ksi2  -1 + ksi1
                 1 - ksi2  -1 - ksi1
                 1 + ksi2   1 + ksi1
                -1 - ksi2   1 - ksi1])


  f = Rad * (XT1^2 + XT2^2 + XT3^2)^(-3/2)
  dX1dXT1 = f * (XT2^2 + XT3^2)
  dX1dXT2 = -f * XT1 * XT2 
  dX1dXT3 = -f * XT1 * XT3
  dX2dXT1 = dX1dXT2 
  dX2dXT2 = f *(XT1^2 + XT3^2)
  dX2dXT3 = -f * XT2 * XT3
  dX3dXT1 = dX1dXT3 
  dX3dXT2 = dX2dXT3 
  dX3dXT3 = f * (XT1^2 + XT2^2)

  J1 = @SArray([dX1dXT1 dX1dXT2 dX1dXT3   
                dX2dXT1 dX2dXT2 dX2dXT3
                dX3dXT1 dX3dXT2 dX3dXT3])   
  J .= 0.25 * J1 * JP * J3      

  X .= XLoc / norm(XLoc) * Rad
end

@inline function Jacobi!(J,detJ,pinvJ,X,::Grids.Quad,ksi1,ksi2,P1,P2,P3,P4,Rad)

  XT1 =  0.25*(P1.x*(1-ksi1)*(1-ksi2)+
            P2.x*(1+ksi1)*(1-ksi2)+
            P3.x*(1+ksi1)*(1+ksi2)+
            P4.x*(1-ksi1)*(1+ksi2))
    
  XT2 =  0.25*(P1.y*(1-ksi1)*(1-ksi2)+
            P2.y*(1+ksi1)*(1-ksi2)+
            P3.y*(1+ksi1)*(1+ksi2)+
            P4.y*(1-ksi1)*(1+ksi2))
           
  XT3 =  0.25*(P1.z*(1-ksi1)*(1-ksi2)+
           P2.z*(1+ksi1)*(1-ksi2)+
            P3.z*(1+ksi1)*(1+ksi2)+
            P4.z*(1-ksi1)*(1+ksi2))
    
  XLoc = SVector{3}(XT1,XT2,XT3)
  JP = @SArray[P1.x P2.x P3.x P4.x;
               P1.y P2.y P3.y P4.y;
               P1.z P2.z P3.z P4.z]

  J3 = @SArray([-0.25 + 0.25*ksi2  -0.25 + 0.25*ksi1
                 0.25 - 0.25*ksi2  -0.25 - 0.25*ksi1
                 0.25 + 0.25*ksi2   0.25 + 0.25*ksi1
                -0.25 - 0.25*ksi2   0.25 - 0.25*ksi1])


  f = Rad *(XT1^2 + XT2^2 + XT3^2)^(-3/2)
  dX1dXT1 = f * (XT2^2 + XT3^2)
  dX1dXT2= -f * XT1 * XT2 
  dX1dXT3= -f * XT1 * XT3
  dX2dXT1 = dX1dXT2 
  dX2dXT2 = f * (XT1^2+XT3^2)
  dX2dXT3 = -f * XT2 * XT3
  dX3dXT1 = dX1dXT3 
  dX3dXT2 = dX2dXT3 
  dX3dXT3 = f*(XT1^2+XT2^2)

  J1  =   @SArray([dX1dXT1    dX1dXT2     dX1dXT3   
                dX2dXT1     dX2dXT2     dX2dXT3
                dX3dXT1     dX3dXT2     dX3dXT3])   
  J   .=   J1*JP*J3

  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ  .= pinvJac(J)
  X .= XLoc / norm(XLoc) * Rad
end

@inline function Jacobi!(J,X,::Grids.Quad,ksi1,ksi2,P1,P2,P3,P4,Rad)

  XT1 =  0.25*(P1.x*(1-ksi1)*(1-ksi2)+
            P2.x*(1+ksi1)*(1-ksi2)+
            P3.x*(1+ksi1)*(1+ksi2)+
            P4.x*(1-ksi1)*(1+ksi2))
    
  XT2 =  0.25*(P1.y*(1-ksi1)*(1-ksi2)+
            P2.y*(1+ksi1)*(1-ksi2)+
            P3.y*(1+ksi1)*(1+ksi2)+
            P4.y*(1-ksi1)*(1+ksi2))
           
  XT3 =  0.25*(P1.z*(1-ksi1)*(1-ksi2)+
           P2.z*(1+ksi1)*(1-ksi2)+
            P3.z*(1+ksi1)*(1+ksi2)+
            P4.z*(1-ksi1)*(1+ksi2))
    
  XLoc = SVector{3}(XT1,XT2,XT3)
  JP = @SArray[P1.x P2.x P3.x P4.x;
               P1.y P2.y P3.y P4.y;
               P1.z P2.z P3.z P4.z]

  J3 = @SArray([-1 + ksi2  -1 + ksi1
                 1 - ksi2  -1 - ksi1
                 1 + ksi2   1 + ksi1
                -1 - ksi2   1 - ksi1])

  f = Rad *(XT1^2 + XT2^2 + XT3^2)^(-3/2)
  dX1dXT1 = f * (XT2^2 + XT3^2)
  dX1dXT2= -f * XT1 * XT2 
  dX1dXT3= -f * XT1 * XT3
  dX2dXT1 = dX1dXT2 
  dX2dXT2 = f * (XT1^2 + XT3^2)
  dX2dXT3 = -f * XT2 * XT3
  dX3dXT1 = dX1dXT3 
  dX3dXT2 = dX2dXT3 
  dX3dXT3 = f*(XT1^2 + XT2^2)

  J1 = @SArray([dX1dXT1 dX1dXT2 dX1dXT3   
                dX2dXT1 dX2dXT2 dX2dXT3
                dX3dXT1 dX3dXT2 dX3dXT3])   
  J .= 0.25 * J1 * JP * J3

  X .= XLoc / norm(XLoc) * Rad
end

function JacobiCart!(J, detJ, pinvJ, X, ::Grids.Tri, ksi1, ksi2, F, Grid)
  XT1 = 0.5 * (F.P[1].x * (-ksi1 - ksi2) +
               F.P[2].x * (1 + ksi1) +
               F.P[3].x * (1 + ksi2))

  XT2 = 0.5 * (F.P[1].y * (-ksi1 - ksi2) +
               F.P[2].y * (1 + ksi1) +
               F.P[3].y * (1 + ksi2))

  XT3 = 0.5 * (F.P[1].z * (-ksi1 - ksi2) +
               F.P[2].z * (1 + ksi1) +
               F.P[3].z * (1 + ksi2))
    
  XLoc = SVector{3}(XT1,XT2,XT3)

  JP = @SArray[F.P[1].x F.P[2].x F.P[3].x;
               F.P[1].y F.P[2].y F.P[3].y;
               F.P[1].z F.P[2].z F.P[3].z]
  
  J3 = @SArray([-1/2 -1/2;
                 1/2   0;
                  0   1/2])
  J .= JP * J3

  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ .= pinvJac(J)
  X .= XLoc
end

function JacobiCart!(J, detJ, pinvJ, X, ::Grids.Quad, ksi1, ksi2, F, Grid)
  XT1 = 0.25 * (F.P[1].x * (1 - ksi1) * (1 - ksi2) +
                F.P[2].x * (1 + ksi1) * (1 - ksi2) +
                F.P[3].x * (1 + ksi1) * (1 + ksi2) +
                F.P[4].x * (1 - ksi1) * (1 + ksi2))

  XT2 = 0.25 * (F.P[1].y * (1 - ksi1) * (1 - ksi2) +
                F.P[2].y * (1 + ksi1) * (1 - ksi2) +
                F.P[3].y * (1 + ksi1) * (1 + ksi2) +
                F.P[4].y * (1 - ksi1) * (1 + ksi2))

  XT3 = 0.25 * (F.P[1].z * (1 - ksi1) * (1 - ksi2) +
                F.P[2].z * (1 + ksi1) * (1 - ksi2) +
                F.P[3].z * (1 + ksi1) * (1 + ksi2) +
                F.P[4].z * (1 - ksi1) * (1 + ksi2))

  XLoc = SVector{3}(XT1,XT2,XT3)

  JP = @SArray[F.P[1].x F.P[2].x F.P[3].x F.P[4].x;
                 F.P[1].y F.P[2].y F.P[3].y F.P[4].y;
                 F.P[1].z F.P[2].z F.P[3].z F.P[4].z]

  J3 = @SArray([-0.25 + 0.25 * ksi2  -0.25 + 0.25 * ksi1
                 0.25 - 0.25 * ksi2  -0.25 - 0.25 * ksi1
                 0.25 + 0.25 * ksi2   0.25 + 0.25 * ksi1
                -0.25 - 0.25 * ksi2   0.25 - 0.25 * ksi1])

  J .= JP * J3

  @views detJLoc = det(J[:,1],J[:,2])
  detJ .= detJLoc
  pinvJ .= pinvJac(J)
  X .= XLoc
end

@inline function det(a,b)

  d = (a[2] * b[3] - a[3] * b[2])^2 +
      (a[1] * b[3] - a[3] * b[1])^2 +
      (a[1] * b[2] - a[2] * b[1])^2 
  d = sqrt(d)
end  

@inline function pinvJac(J)
  g11 = J[1,1] * J[1,1] + J[2,1] * J[2,1] + J[3,1] * J[3,1]
  g12 = J[1,1] * J[1,2] + J[2,1] * J[2,2] + J[3,1] * J[3,2]
  g22 = J[1,2] * J[1,2] + J[2,2] * J[2,2] + J[3,2] * J[3,2]
  det = g11 * g22 - g12^2
  i11 = g22 / det
  i21 = -g12 / det
  i12 = -g12 / det
  i22 = g11 / det
  pJ11 = J[1,1] * i11 + J[1,2] * i21
  pJ12 = J[1,1] * i12 + J[1,2] * i22
  pJ21 = J[2,1] * i11 + J[2,2] * i21
  pJ22 = J[2,1] * i12 + J[2,2] * i22
  pJ31 = J[3,1] * i11 + J[3,2] * i21
  pJ32 = J[3,1] * i12 + J[3,2] * i22
  @SArray([pJ11 pJ12
          pJ21 pJ22
          pJ31 pJ32])
   
end
