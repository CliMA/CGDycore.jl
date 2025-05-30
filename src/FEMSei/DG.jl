function Lagrange(x,xw,i)
  L = 1 + 0*x
  @inbounds for j = 1 : length(xw)
    if j != i
      L = L * (x - xw[j]) / (xw[i] - xw[j])
    end
  end
  return L
end

mutable struct DGStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: ScalarElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2} 
  Gradphi::Array{Polynomial,3}                           
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  points::Array{Float64,2}
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

#DG Quad

function DGStruct{FT}(backend,k::Int,Type::Grids.Quad,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  kp1 = k + 1
  DoF = kp1 * kp1
  Comp = 1
  points = KernelAbstractions.zeros(backend,Float64,DoF,2)
  phi = Array{Polynomial,2}(undef,DoF,1)
  Gradphi = Array{Polynomial,3}(undef,DoF,1,2)
  L1 = Array{Polynomial,1}(undef,kp1)
  L2 = Array{Polynomial,1}(undef,kp1)

  @polyvar x[1:2] 
  if k == 0
    phi[1,1] = 1.0 + 0.0*x[1] + 0.0*x[2]
    xw = zeros(1)
  else
    xw,_= gausslobatto(kp1)
  end
  @inbounds for i = 1 : kp1
    L1[i] = Lagrange(x[1],xw,i)   
    L2[i] = Lagrange(x[2],xw,i)  
  end  
  iDoF = 1
  @inbounds for j = 1 : kp1
    @inbounds for i = 1 : kp1
      phi[iDoF,1] = L1[i] * L2[j]
      points[iDoF,1] = xw[i]
      points[iDoF,2] = xw[j]
      iDoF += 1
    end  
  end  
  @inbounds for iDoF = 1 : DoF
    Gradphi[iDoF,1,1] = differentiate(phi[iDoF,1],x[1])
    Gradphi[iDoF,1,2] = differentiate(phi[iDoF,1],x[2])
  end
  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = DoF * Grid.NumFaces
  NumI = DoF * Grid.NumFaces
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : DoF
      GlobCPU[iD,iF] = DoF * (Grid.Faces[iF].F - 1) + iD
    end
  end
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return DGStruct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Gradphi,
    NumG,
    NumI,
    Type,
    points,
    M,
    LUM,
      )
end

#DG Tri

function DGStruct{FT}(backend,k,Type::Grids.Tri,Grid) where FT<:AbstractFloat

  if k == 0
    DoF = 1
    Comp = 1
    @polyvar x1 x2
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    phi[1,1] = 1.0 + 0.0*x1 + 0.0*x2
    @inbounds for i = 1 : DoF
      @inbounds for j = 1 : Comp
        Gradphi[i,j,1] = differentiate(phi[i,j],x1)
        Gradphi[i,j,2] = differentiate(phi[i,j],x2)
      end
    end 
    points = KernelAbstractions.zeros(backend,Float64,DoF,2)
    points = [-1/3 -1/3]     
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = Grid.NumFaces
    NumI = Grid.NumFaces
    @inbounds for iF = 1 : Grid.NumFaces
      GlobCPU[1,iF] = Grid.Faces[iF].F
    end
  else
    Comp = 1  
    DoF, DoFE, DoFF, phi, Gradphi, points = FEMSei.ConstructCG(k,Type)  
    NumG = DoF * Grid.NumFaces 
    NumI = NumG
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    @inbounds for iF = 1 : Grid.NumFaces
      @inbounds for iDoF = 1 : DoF
        GlobCPU[iDoF,iF] = iDoF + (Grid.Faces[iF].F - 1) * DoF
      end  
    end
  end  
  Glob = KernelAbstractions.zeros(backend,Int,size(GlobCPU))
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)

  return DGStruct{FT,
                  typeof(Glob)}( 
    Glob,
    DoF,
    Comp,
    phi,
    Gradphi,
    NumG,
    NumI,
    Type,
    points,
    M,
    LUM,
      )
end
