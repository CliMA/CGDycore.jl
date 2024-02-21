  import CGDycore:
  Examples, Parallels, Models, Grids, Outputs, Integration,  GPU, DyCore, FEMSei, FiniteVolumes
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using DynamicPolynomials
using FastGaussQuadrature

FTB = Float32
QuadOrd = 2
backend = CPU()

QQ = FEMSei.QuadRule{FTB}(Grids.Quad(),backend,QuadOrd)

function IntegrateFace(phi,psi,QQ)
  Weights = QQ.Weights
  Points = QQ.Points
  Int = 0.0
  for i = 1 : length(Weights)
     Int += Weights[i]*phi(Points[i,1],Points[i,2])*psi(Points[i,1],Points[i,2]) 
  end   
  return Int
end

function IntegrateEdgeX(phi,psi,X)
  x,w = gausslobatto(3)
  Int = 0.0
  for i = 1 : length(w)
    Int += w[i] * phi(X,x[i]) * psi(X,x[i])     
  end
  return Int
end  

function IntegrateEdgeY(phi,psi,Y)
  x,w = gausslobatto(3)
  Int = 0.0
  for i = 1 : length(w)
    Int += w[i] * phi(x[i],Y) * psi(x[i],Y)     
  end
  return Int
end  


  DoF = 4
  Comp = 1
  @polyvar x y
  phi = Array{Polynomial,2}(undef,DoF,Comp)
  xP, w = gaussradau(2)
  lx0 = (x - xP[2])/(xP[1] - xP[2])
  lx1 = (x - xP[1])/(xP[2] - xP[1])
  ly0 = (y - xP[2])/(xP[1] - xP[2])
  ly1 = (y - xP[1])/(xP[2] - xP[1])
  phi[1,1] = lx0 * ly0
  phi[2,1] = lx1 * ly0
  phi[3,1] = lx0 * ly1
  phi[4,1] = lx1 * ly1
  Gradphi = Array{Polynomial,2}(undef,DoF,2)
  for i = 1 : DoF
    Gradphi[i,1] = differentiate(phi[i,1],x)
    Gradphi[i,2] = differentiate(phi[i,1],y)
  end


  DoF = 8
  Comp = 2
  @polyvar x y
  psi = Array{Polynomial,2}(undef,DoF,Comp)
  Divpsi = Array{Polynomial,2}(undef,DoF,1)
  xP, w = gaussradau(2)
  xP .= -xP
  lx0 = (x - xP[2])/(xP[1] - xP[2])
  lx1 = (x - xP[1])/(xP[2] - xP[1])
  ly0 = (y - xP[2])/(xP[1] - xP[2])
  ly1 = (y - xP[1])/(xP[2] - xP[1])
  p0 = 0.0*x + 0.0*y
  psi[1,1] = p0
  psi[1,2] = lx0 * ly0
  psi[2,1] = p0
  psi[2,2] = lx1 * ly0

  psi[3,1] = lx0 * ly0
  psi[3,2] = p0
  psi[4,1] = lx0 * ly1
  psi[4,2] = p0

  psi[5,1] = lx1 * ly0
  psi[5,2] = p0
  psi[6,1] = lx1 * ly1
  psi[6,2] = p0

  psi[7,1] = p0
  psi[7,2] = lx0 * ly1
  psi[8,1] = p0
  psi[8,2] = lx1 * ly1

  Divpsi[1,1] =  (differentiate(psi[1,1],x) + differentiate(psi[1,2],y))
  Divpsi[2,1] =  (differentiate(psi[2,1],x) + differentiate(psi[2,2],y))
  Divpsi[3,1] =  (differentiate(psi[3,1],x) + differentiate(psi[3,2],y))
  Divpsi[4,1] =  (differentiate(psi[4,1],x) + differentiate(psi[4,2],y))
  Divpsi[5,1] =  (differentiate(psi[5,1],x) + differentiate(psi[5,2],y))
  Divpsi[6,1] =  (differentiate(psi[6,1],x) + differentiate(psi[6,2],y))
  Divpsi[7,1] =  (differentiate(psi[7,1],x) + differentiate(psi[7,2],y))
  Divpsi[8,1] =  (differentiate(psi[8,1],x) + differentiate(psi[8,2],y))

  Int =zeros(8,4)
# psi[5,1] = lx1 * ly0
# psi[5,2] = p0
  for i = 1 : 8
    for j = 1 : 4  
      Int[i,j] = IntegrateFace(Divpsi[i,1],phi[j,1],QQ) +
        IntegrateEdgeX(psi[i,1],phi[j,1],-1.0) + 
        IntegrateEdgeY(psi[i,2],phi[j,1],-1.0) 
    end    
  end  
#=
  Int[5,2] = IntegrateFace(Divpsi[5,1],phi[2,1],QQ) +
    IntegrateEdgeX(psi[5,1],phi[2,1],-1.0) + 
    IntegrateEdgeY(psi[5,2],phi[2,1],-1.0) 
  Int[5,3] = IntegrateFace(Divpsi[5,1],phi[3,1],QQ) +
    IntegrateEdgeX(psi[5,1],phi[3,1],-1.0) + 
    IntegrateEdgeY(psi[5,2],phi[3,1],-1.0) 
  Int[5,4] = IntegrateFace(Divpsi[5,1],phi[4,1],QQ) +
    IntegrateEdgeX(psi[5,1],phi[4,1],-1.0) + 
    IntegrateEdgeY(psi[5,2],phi[4,1],-1.0) 
  Int[6,1] = IntegrateFace(Divpsi[6,1],phi[1,1],QQ) +
    IntegrateEdgeX(psi[6,1],phi[1,1],-1.0) + 
    IntegrateEdgeY(psi[6,2],phi[1,1],-1.0) 
  Int[6,2] = IntegrateFace(Divpsi[6,1],phi[2,1],QQ) +
    IntegrateEdgeX(psi[6,1],phi[2,1],-1.0) + 
    IntegrateEdgeY(psi[6,2],phi[2,1],-1.0) 
  Int[6,3] = IntegrateFace(Divpsi[6,1],phi[3,1],QQ) +
    IntegrateEdgeX(psi[6,1],phi[3,1],-1.0) + 
    IntegrateEdgeY(psi[6,2],phi[3,1],-1.0) 
  Int[6,4] = IntegrateFace(Divpsi[6,1],phi[4,1],QQ) +
    IntegrateEdgeX(psi[6,1],phi[4,1],-1.0) + 
    IntegrateEdgeY(psi[6,2],phi[4,1],-1.0) 
# psi[7,1] = p0
# psi[7,2] = lx0 * ly1  
  Int[7,1] = IntegrateFace(Divpsi[7,1],phi[1,1],QQ) +
    IntegrateEdgeX(psi[7,1],phi[1,1],-1.0) + 
    IntegrateEdgeY(psi[7,2],phi[1,1],-1.0)
  Int[7,2] = IntegrateFace(Divpsi[7,1],phi[2,1],QQ) +
    IntegrateEdgeX(psi[7,1],phi[2,1],-1.0) + 
    IntegrateEdgeY(psi[7,2],phi[2,1],-1.0)
  Int[7,3] = IntegrateFace(Divpsi[7,1],phi[3,1],QQ) +
    IntegrateEdgeX(psi[7,1],phi[3,1],-1.0) + 
    IntegrateEdgeY(psi[7,2],phi[3,1],-1.0)
  Int[7,4] = IntegrateFace(Divpsi[7,1],phi[4,1],QQ) +
    IntegrateEdgeX(psi[7,1],phi[4,1],-1.0) + 
    IntegrateEdgeY(psi[7,2],phi[4,1],-1.0)

  Int[8,1] = IntegrateFace(Divpsi[8,1],phi[1,1],QQ) +
    IntegrateEdgeX(psi[8,1],phi[1,1],-1.0) + 
    IntegrateEdgeY(psi[8,2],phi[1,1],-1.0)
  Int[8,2] = IntegrateFace(Divpsi[8,1],phi[2,1],QQ) +
    IntegrateEdgeX(psi[8,1],phi[2,1],-1.0) + 
    IntegrateEdgeY(psi[8,2],phi[2,1],-1.0)
  Int[8,3] = IntegrateFace(Divpsi[8,1],phi[3,1],QQ) +
    IntegrateEdgeX(psi[8,1],phi[3,1],-1.0) + 
    IntegrateEdgeY(psi[8,2],phi[3,1],-1.0)
  Int[8,4] = IntegrateFace(Divpsi[8,1],phi[4,1],QQ) +
    IntegrateEdgeX(psi[8,1],phi[4,1],-1.0) + 
    IntegrateEdgeY(psi[8,2],phi[4,1],-1.0)
=#   
