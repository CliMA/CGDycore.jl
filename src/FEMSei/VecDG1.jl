mutable struct VecDGStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: VectorElement
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

#VecDG Quad
function VecDGStruct{FT}(backend,k::Int,Type::Grids.ElementType,Grid) where FT<:AbstractFloat
    
    if k == 0 #order 0
        if Type == Grids.Tri() #Tri
            Glob = KernelAbstractions.zeros(backend,Int,0,0)
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
            points = KernelAbstractions.zeros(backend,Float64,DoF,2)
            points = [0.0 0.0]
          
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

        else Type == Grids.Quad() #Quad
            Glob = KernelAbstractions.zeros(backend,Int,0,0)
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
            points = KernelAbstractions.zeros(backend,Float64,DoF,2)
            points = [0.0 0.0]
      
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
        end

    elseif k == 1 #order 1
        if Type == Grids.Tri() #Tri
           points = KernelAbstractions.zeros(backend,Float64,DoF,2)
            Glob = KernelAbstractions.zeros(backend,Int,0,0)
            DoF = 9 
            Comp = 3
            @polyvar x1 x2
            phi = Array{Polynomial,2}(undef,DoF,Comp)
            phi[1,1] = 0.0*x1 + 0.0*x2 + 2.0
            phi[1,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[1,3] = 0.0 + 0.0*x1 + 0.0*x2
    
            phi[2,1] = 0.0*x1 + 0.0*x2 + 2.0
            phi[2,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[2,3] = 0.0 + 0.0*x1 + 0.0*x2
    
            phi[3,1] = 0.0*x1 + 0.0*x2 + 2.0
            phi[3,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[3,3] = 0.0 + 0.0*x1 + 0.0*x2
    
            phi[4,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[4,2] = 0.0*x1 + 1.0*x2 - 1.0/3.0
            phi[4,3] = 0.0 + 0.0*x1 + 0.0*x2
    
            phi[5,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[5,2] = 0.0*x1 + 1.0*x2 - 1.0/3.0
            phi[5,3] = 0.0 + 0.0*x1 + 0.0*x2
    
            phi[6,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[6,2] = 0.0*x1 + 1.0*x2 - 1.0/3.0
            phi[6,3] = 0.0 + 0.0*x1 + 0.0*x2
    
            phi[7,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[7,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[7,3] = 1.0*x1 + 0.0*x2 - 1.0/3.0
    
            phi[8,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[8,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[8,3] = 1.0*x1 + 0.0*x2 - 1.0/3.0
    
            phi[9,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[9,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[9,3] = 1.0*x1 + 0.0*x2 - 1.0/3.0
    
           
            Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
            for i = 1 : DoF
                for j = 1 : Comp
                    Gradphi[i,j,1] = differentiate(phi[i,j],x1)
                    Gradphi[i,j,2] = differentiate(phi[i,j],x2)
                end
            end
    
            points = KernelAbstractions.zeros(backend,Float64,3,2)
            points[1,:] = [-1.0 -1.0]
            points[2,:] = [1.0 -1.0]
            points[3,:] = [-1.0 1.0]
    
            Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
            @show size(Glob)
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

        else Type == Grids.Quad() #Quad
            Glob = KernelAbstractions.zeros(backend,Int,0,0)
            DoF = 12
            Comp = 3
            @polyvar x1 x2
            phi = Array{Polynomial,2}(undef,DoF,Comp)
            phi[1,1] = 1/4 * (1.0-1.0*x1) * (1.0-1.0*x2)
            phi[1,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[1,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[2,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[2,2] = 1/4 * (1.0-1.0*x1) * (1.0-1.0*x2)
            phi[2,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[3,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[3,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[3,3] = 1/4 * (1.0-1.0*x1) * (1.0-1.0*x2)

            phi[4,1] = 1/4 * (1.0+1.0*x1) * (1.0-1.0*x2)
            phi[4,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[4,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[5,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[5,2] = 1/4 * (1.0+1.0*x1) * (1.0-1.0*x2)
            phi[5,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[6,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[6,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[6,3] = 1/4 * (1.0+1.0*x1) * (1.0-1.0*x2)

            phi[7,1] = 1/4 * (1.0+1.0*x1) * (1.0+1.0*x2)
            phi[7,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[7,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[8,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[8,2] = 1/4 * (1.0+1.0*x1) * (1.0+1.0*x2)
            phi[8,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[9,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[9,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[9,3] = 1/4 * (1.0+1.0*x1) * (1.0+1.0*x2)

            phi[10,1] = 1/4 * (1.0-1.0*x1) * (1.0+1.0*x2)
            phi[10,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[10,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[11,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[11,2] = 1/4 * (1.0-1.0*x1) * (1.0+1.0*x2)
            phi[11,3] = 0.0 + 0.0*x1 + 0.0*x2

            phi[12,1] = 0.0 + 0.0*x1 + 0.0*x2
            phi[12,2] = 0.0 + 0.0*x1 + 0.0*x2
            phi[12,3] = 1/4 * (1.0-1.0*x1) * (1.0+1.0*x2)

            Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
            for i = 1 : DoF
                for j = 1 : Comp
                    Gradphi[i,j,1] = differentiate(phi[i,j],x1)
                    Gradphi[i,j,2] = differentiate(phi[i,j],x2)
                end
            end

            points = KernelAbstractions.zeros(backend,Float64,4,2)
            points[1,:] = [-1.0 -1.0]
            points[2,:] = [1.0 -1.0]
            points[3,:] = [1.0 1.0]
            points[4,:] = [-1.0 1.0]
      
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
        end

    else println("not defined")
    end

    return VecDGStruct{FT,
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