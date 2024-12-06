mutable struct VecDG0Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: VectorElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2} 
  Gradphi::Array{Polynomial,3}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

#VecDG0 Quad

function VecDG0Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
    Glob = KernelAbstractions.zeros(backend,Int,0,0)
    Type = Grids.Quad()
    DoF = 3 #vorn
    Comp = 3
    @polyvar x1 x2
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    phi[1,1] = 1.0 + 0.0*x1 + 0.0*x2
    phi[1,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[1,3] = 0.0 + 0.0*x1 + 0.0*x2
    phi[2,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[2,2] = 1.0 + 0.0*x1 + 0.0*x2
    phi[2,3] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,3] = 1.0 + 0.0*x1 + 0.0*x2
    
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    for i = 1 : DoF
        for j = 1 : Comp
            Gradphi[i,j,1] = differentiate(phi[i,j],x1)
            Gradphi[i,j,2] = differentiate(phi[i,j],x2)
        end
    end
      
    Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = DoF * Grid.NumFaces
    NumI = DoF * Grid.NumFaces
    for iF = 1 : Grid.NumFaces
        for iD = 1 : DoF
            GlobCPU[iD,iF] = (Grid.Faces[iF].F - 1 ) * DoF + iD
        end
    end
    copyto!(Glob,GlobCPU)
    M = sparse([1],[1],[1.0])
    LUM = lu(M)
    return VecDG0Struct{FT,
                    typeof(Glob)}( 
      Glob,
      DoF,
      Comp,
      phi,
      Gradphi,                      
      NumG,
      NumI,
      Type,
      M,
      LUM,
        )
end

#VecDG0 Tri

function VecDG0Struct{FT}(::Grids.Tri,backend,Grid) where FT<:AbstractFloat
    Glob = KernelAbstractions.zeros(backend,Int,0,0)
    Type = Grids.Tri()
    DoF = 3 
    Comp = 3
    @polyvar x1 x2
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    phi[1,1] = 1.0 + 0.0*x1 + 0.0*x2
    phi[1,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[1,3] = 0.0 + 0.0*x1 + 0.0*x2
    phi[2,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[2,2] = 1.0 + 0.0*x1 + 0.0*x2
    phi[2,3] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,3] = 1.0 + 0.0*x1 + 0.0*x2
    
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    for i = 1 : DoF
        for j = 1 : Comp
            Gradphi[i,j,1] = differentiate(phi[i,j],x1)
            Gradphi[i,j,2] = differentiate(phi[i,j],x2)
        end
    end
      
    Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = DoF * Grid.NumFaces
    NumI = DoF * Grid.NumFaces
    for iF = 1 : Grid.NumFaces
        for iD = 1 : DoF
            GlobCPU[iD,iF] = (Grid.Faces[iF].F - 1 ) * DoF + iD
        end
    end
    copyto!(Glob,GlobCPU)
    M = sparse([1],[1],[1.0])
    LUM = lu(M)
    return VecDG0Struct{FT,
                    typeof(Glob)}( 
      Glob,
      DoF,
      Comp,
      phi,
      Gradphi,                      
      NumG,
      NumI,
      Type,
      M,
      LUM,
        )
end

mutable struct VecDG1Struct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: VectorElement
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2} 
  Gradphi::Array{Polynomial,3}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

#VecDG1 Quad

function VecDG1Struct{FT}(::Grids.Quad,backend,Grid) where FT<:AbstractFloat
    Glob = KernelAbstractions.zeros(backend,Int,0,0)
    Type = Grids.Quad()
    DoF = 12 #vorn
    Comp = 3
    @polyvar x1 x2
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    phi[1,1] = (1.0-1.0*x1) * (1.0-1.0*x2)
    phi[1,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[1,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[2,1] = (1.0+1.0*x1) * (1.0-1.0*x2)
    phi[2,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[2,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[3,1] = (1.0+1.0*x1) * (1.0+1.0*x2)
    phi[3,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[3,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[4,1] = (1.0-1.0*x1) * (1.0+1.0*x2)
    phi[4,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[4,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[5,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[5,2] = (1.0-1.0*x1) * (1.0-1.0*x2)
    phi[5,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[6,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[6,2] = (1.0+1.0*x1) * (1.0-1.0*x2)
    phi[6,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[7,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[7,2] = (1.0+1.0*x1) * (1.0+1.0*x2)
    phi[7,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[8,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[8,2] = (1.0-1.0*x1) * (1.0+1.0*x2)
    phi[8,3] = 0.0 + 0.0*x1 + 0.0*x2

    phi[9,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[9,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[9,3] = (1.0-1.0*x1) * (1.0-1.0*x2)

    phi[10,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[10,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[10,3] = (1.0+1.0*x1) * (1.0-1.0*x2)

    phi[11,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[11,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[11,3] = (1.0+1.0*x1) * (1.0+1.0*x2)

    phi[12,1] = 0.0 + 0.0*x1 + 0.0*x2
    phi[12,2] = 0.0 + 0.0*x1 + 0.0*x2
    phi[12,3] = (1.0-1.0*x1) * (1.0+1.0*x2)
    
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    for i = 1 : DoF
        for j = 1 : Comp
            Gradphi[i,j,1] = differentiate(phi[i,j],x1)
            Gradphi[i,j,2] = differentiate(phi[i,j],x2)
        end
    end
      
    Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = DoF * Grid.NumFaces
    NumI = DoF * Grid.NumFaces
    for iF = 1 : Grid.NumFaces
        for iD = 1 : DoF
            GlobCPU[iD,iF] = (Grid.Faces[iF].F - 1 ) * DoF + iD
        end
    end
    copyto!(Glob,GlobCPU)
    M = sparse([1],[1],[1.0])
    LUM = lu(M)
    return VecDG1Struct{FT,
                    typeof(Glob)}( 
      Glob,
      DoF,
      Comp,
      phi,
      Gradphi,                      
      NumG,
      NumI,
      Type,
      M,
      LUM,
        )
end

#VecDG1 Tri

function VecDG1Struct{FT}(::Grids.Tri,backend,Grid) where FT<:AbstractFloat
    Glob = KernelAbstractions.zeros(backend,Int,0,0)
    Type = Grids.Tri()
    DoF = 9 #vorn
    Comp = 3
    @polyvar ksi1 ksi2
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    phi[1,1] = 0.0*ksi1 + 0.0*ksi2 + 2.0
    phi[1,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[1,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2

    phi[2,1] = 0.0*ksi1 + 0.0*ksi2 + 2.0
    phi[2,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[2,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2

    phi[3,1] = 0.0*ksi1 + 0.0*ksi2 + 2.0
    phi[3,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[3,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2

    phi[4,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[4,2] = 0.0*ksi1 + 1.0*ksi2 - 1.0/3.0
    phi[4,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2

    phi[5,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[5,2] = 0.0*ksi1 + 1.0*ksi2 - 1.0/3.0
    phi[5,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2

    phi[6,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[6,2] = 0.0*ksi1 + 1.0*ksi2 - 1.0/3.0
    phi[6,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2

    phi[7,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[7,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[7,3] = 1.0*ksi1 + 0.0*ksi2 - 1.0/3.0

    phi[8,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[8,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[8,3] = 1.0*ksi1 + 0.0*ksi2 - 1.0/3.0

    phi[9,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[9,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    phi[9,3] = 1.0*ksi1 + 0.0*ksi2 - 1.0/3.0

    for s = 1 : DoF
        for t = 1 : Comp
          phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
        end
    end
    
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    for i = 1 : DoF
        for j = 1 : Comp
            Gradphi[i,j,1] = differentiate(phi[i,j],x1)
            Gradphi[i,j,2] = differentiate(phi[i,j],x2)
        end
    end
      
    Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
    GlobCPU = zeros(Int,DoF,Grid.NumFaces)
    NumG = DoF * Grid.NumFaces
    NumI = DoF * Grid.NumFaces
    for iF = 1 : Grid.NumFaces
        for iD = 1 : DoF
            GlobCPU[iD,iF] = (Grid.Faces[iF].F - 1 ) * DoF + iD
        end
    end
    copyto!(Glob,GlobCPU)
    M = sparse([1],[1],[1.0])
    LUM = lu(M)
    return VecDG1Struct{FT,
                    typeof(Glob)}( 
      Glob,
      DoF,
      Comp,
      phi,
      Gradphi,                      
      NumG,
      NumI,
      Type,
      M,
      LUM,
        )
end