mutable struct CGStruct{FT<:AbstractFloat,
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

#CG Quad

function CGStruct{FT}(backend,k::Int,Type::Grids.Quad,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  kp1 = k + 1
  km1 = k - 1
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
  NumG = (k - 1) * (k - 1) * Grid.NumFaces +
    (k - 1) * Grid.NumEdges + Grid.NumNodes
  NumI = NumG
# Ordnung Ansatzfunction Freiheitsgrade
  @inbounds for iF = 1 : Grid.NumFaces
    iD = 1
    GlobCPU[iD,iF] = Grid.Faces[iF].N[1]
    iD = kp1
    GlobCPU[iD,iF] = Grid.Faces[iF].N[2]
    iD = k * kp1 + 1
    GlobCPU[iD,iF] = Grid.Faces[iF].N[4]
    iD = kp1 * kp1
    GlobCPU[iD,iF] = Grid.Faces[iF].N[3]
    for i = 1 : km1
      iD = 1 + i  
      GlobCPU[iD,iF] = i + (Grid.Faces[iF].E[1] - 1) * km1 + Grid.NumNodes
      iD = 1 + i + k * kp1  
      GlobCPU[iD,iF] = i + (Grid.Faces[iF].E[3] - 1) * km1 + Grid.NumNodes
      iD = 1 + i * kp1  
      GlobCPU[iD,iF] = i + (Grid.Faces[iF].E[4] - 1) * km1 + Grid.NumNodes
      iD = (i + 1) * kp1  
      GlobCPU[iD,iF] = i + (Grid.Faces[iF].E[2] - 1) * km1 + Grid.NumNodes
    end  
    @inbounds for j = 1 : km1
      @inbounds for i = 1 : km1
        iD = j * kp1 + i + 1
        GlobCPU[iD,iF] = i + (j - 1) * km1 + km1 * km1 * (Grid.Faces[iF].F - 1) +
          Grid.NumEdges * km1 + Grid.NumNodes
      end  
    end
  end
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  return CGStruct{FT,
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

#CG Tri

function CGStruct{FT}(backend,k,Type::Grids.Tri,Grid) where FT<:AbstractFloat
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
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
    Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = Grid.NumFaces
    NumI = Grid.NumFaces
    @inbounds for iF = 1 : Grid.NumFaces
      GlobCPU[1,iF] = Grid.Faces[iF].F
    end
    copyto!(Glob,GlobCPU)
    M = sparse([1],[1],[1.0])
    LUM = lu(M)
  elseif k == 1
    DoF = 3
    Comp = 1
    nu = Array{Polynomial,2}(undef,DoF,Comp)
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    @polyvar x1 x2 ksi1 ksi2
    nu[1,1] = -1.0*ksi1 + -1.0*ksi2 + 1.0
    nu[2,1] = 1.0*ksi1 + 0.0*ksi2 + 0.0
    nu[3,1] = 0.0*ksi1 + 1.0*ksi2 + 0.0
  
    @inbounds for s = 1 : DoF
      @inbounds for t = 1 : 1
        phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
      end
    end
    @inbounds for i = 1 : DoF
        @inbounds for j = 1 : Comp
            Gradphi[i,j,1] = differentiate(phi[i,j],x1)
            Gradphi[i,j,2] = differentiate(phi[i,j],x2)
        end
    end
    points = KernelAbstractions.zeros(backend,Float64,DoF,2)
    points[1,:] = [-1.0 -1.0]
    points[2,:] = [1.0 -1.0]
    points[3,:] = [-1.0 1.0]

    Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = 3 * Grid.NumFaces
    NumI = 3 * Grid.NumFaces
    @inbounds for iF = 1 : Grid.NumFaces
      GlobCPU[1,iF] = 3 * Grid.Faces[iF].F - 2
      GlobCPU[2,iF] = 3 * Grid.Faces[iF].F - 1
      GlobCPU[3,iF] = 3 * Grid.Faces[iF].F
    end
    copyto!(Glob,GlobCPU)
    M = sparse([1],[1],[1.0])
    LUM = lu(M)
  else println("Not implemented")
  end

  return CGStruct{FT,
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
