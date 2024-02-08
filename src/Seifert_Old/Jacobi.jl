#TriSphere   
function Jacobi(::Grids.Tri,ksi1,ksi2,F,Grid)
Rad =   Grid.Rad

XT1  =   0.5*(Grid.Nodes[F.N[1]].P.x*(-ksi1)*(-ksi2)+
        Grid.Nodes[F.N[2]].P.x*(1+ksi1)+
        Grid.Nodes[F.N[3]].P.x*(1+ksi2))

XT2  =   0.5*(Grid.Nodes[F.N[1]].P.y*(-ksi1)*(-ksi2)+
        Grid.Nodes[F.N[2]].P.y*(1+ksi1)+
        Grid.Nodes[F.N[3]].P.y*(1+ksi2))
       
XT3  =   0.5*(Grid.Nodes[F.N[1]].P.z*(-ksi1)*(-ksi2)+
        Grid.Nodes[F.N[2]].P.z*(1+ksi1)+
        Grid.Nodes[F.N[3]].P.z*(1+ksi2))

X   =   SVector{3}(XT1,XT2,XT3)

JP  =   @SArray[Grid.Nodes[F.N[1]].P.x Grid.Nodes[F.N[2]].P.x Grid.Nodes[F.N[3]].P.x;
                Grid.Nodes[F.N[1]].P.y Grid.Nodes[F.N[2]].P.y Grid.Nodes[F.N[3]].P.y;
                Grid.Nodes[F.N[1]].P.z Grid.Nodes[F.N[2]].P.z Grid.Nodes[F.N[3]].P.z]

J3  =   @SArray([-1/2 -1/2;
                    1/2 0;
                    0 1/2])


f=Rad*(XT1^2+XT2^2+XT3^2)^(-3/2)
dX1dXT1=f*(XT2^2+XT3^2)
dX1dXT2=-f*XT1*XT2 
dX1dXT3=-f*XT1*XT3
dX2dXT1=dX1dXT2 
dX2dXT2=f*(XT1^2+XT3^2)
dX2dXT3=-f*XT2*XT3
dX3dXT1=dX1dXT3 
dX3dXT2=dX2dXT3 
dX3dXT3=f*(XT1^2+XT2^2)

J1  =   @SArray([dX1dXT1    dX2dXT1     dX3dXT1   
                dX1dXT2     dX2dXT2     dX3dXT2
                dX1dXT3     dX2dXT3     dX3dXT3])   
J   =   J1*JP*J3      

#=
C=0.5*[-1+ksi2  -1+ksi1
  1-ksi2  -1-ksi1
  1+ksi2   1+ksi1]
=#
detJ=norm(cross(J[:,1],J[:,2]))
return J,detJ,X       
end

#Quad   
function Jacobi(::Grids.Quad,ksi1,ksi2,F,Grid)
    Rad =   Grid.Rad
    
    XT1  =   0.25*(Grid.Nodes[F.N[1]].P.x*(1-ksi1)*(1-ksi2)+
            Grid.Nodes[F.N[2]].P.x*(1+ksi1)*(1-ksi2)+
            Grid.Nodes[F.N[3]].P.x*(1+ksi1)*(1+ksi2)+
            Grid.Nodes[F.N[4]].P.x*(1-ksi1)*(1+ksi2))
    
    XT2  =   0.25*(Grid.Nodes[F.N[1]].P.y*(1-ksi1)*(1-ksi2)+
            Grid.Nodes[F.N[2]].P.y*(1+ksi1)*(1-ksi2)+
            Grid.Nodes[F.N[3]].P.y*(1+ksi1)*(1+ksi2)+
            Grid.Nodes[F.N[4]].P.y*(1-ksi1)*(1+ksi2))
           
    XT3  =   0.25*(Grid.Nodes[F.N[1]].P.z*(1-ksi1)*(1-ksi2)+
            Grid.Nodes[F.N[2]].P.z*(1+ksi1)*(1-ksi2)+
            Grid.Nodes[F.N[3]].P.z*(1+ksi1)*(1+ksi2)+
            Grid.Nodes[F.N[4]].P.z*(1-ksi1)*(1+ksi2))
    
    X   =   SVector{3}(XT1,XT2,XT3)
    JP  =   @SArray[Grid.Nodes[F.N[1]].P.x Grid.Nodes[F.N[2]].P.x Grid.Nodes[F.N[3]].P.x Grid.Nodes[F.N[4]].P.x;
                Grid.Nodes[F.N[1]].P.y Grid.Nodes[F.N[2]].P.y Grid.Nodes[F.N[3]].P.y Grid.Nodes[F.N[4]].P.y;
                Grid.Nodes[F.N[1]].P.z Grid.Nodes[F.N[2]].P.z Grid.Nodes[F.N[3]].P.z Grid.Nodes[F.N[4]].P.z]

    J3  =   @SArray([-1+ksi2  -1+ksi1
                    1-ksi2  -1-ksi1
                     1+ksi2   1+ksi1
                    -1-ksi2   1-ksi1])


f=Rad*(XT1^2+XT2^2+XT3^2)^(-3/2)
dX1dXT1=f*(XT2^2+XT3^2)
dX1dXT2=-f*XT1*XT2 
dX1dXT3=-f*XT1*XT3
dX2dXT1=dX1dXT2 
dX2dXT2=f*(XT1^2+XT3^2)
dX2dXT3=-f*XT2*XT3
dX3dXT1=dX1dXT3 
dX3dXT2=dX2dXT3 
dX3dXT3=f*(XT1^2+XT2^2)

J1  =   @SArray([dX1dXT1    dX2dXT1     dX3dXT1   
                dX1dXT2     dX2dXT2     dX3dXT2
                dX1dXT3     dX2dXT3     dX3dXT3])   
J   =   J1*JP*J3      

detJ=norm(cross(J[:,1],J[:,2]))
return J,detJ,X       
end