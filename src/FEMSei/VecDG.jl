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
        end

    elseif k == 1 #order 1
        if Type == Grids.Tri() #Tri
            DoF = 9 
            Comp = 3
            points = KernelAbstractions.zeros(backend,Float64,3,2)
            Glob = KernelAbstractions.zeros(backend,Int,0,0)

            nu = Array{Polynomial,2}(undef,DoF,Comp)
            phi = Array{Polynomial,2}(undef,DoF,Comp)
            Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
            @polyvar ksi1 ksi2 x1 x2

            nu[1,1] = -1.0*ksi1 + -1.0*ksi2 + 1.0
            nu[1,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[1,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    
            nu[2,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[2,2] = -1.0*ksi1 + -1.0*ksi2 + 1.0
            nu[2,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    
            nu[3,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[3,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[3,3] = -1.0*ksi1 + -1.0*ksi2 + 1.0
    
            nu[4,1] = 1.0*ksi1 + 0.0*ksi2 + 0.0
            nu[4,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[4,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    
            nu[5,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[5,2] = 1.0*ksi1 + 0.0*ksi2 + 0.0
            nu[5,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    
            nu[6,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[6,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[6,3] = 1.0*ksi1 + 0.0*ksi2 + 0.0
    
            nu[7,1] = 0.0*ksi1 + 1.0*ksi2 + 0.0
            nu[7,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[7,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    
            nu[8,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[8,2] = 0.0*ksi1 + 1.0*ksi2 + 0.0
            nu[8,3] = 0.0 + 0.0*ksi1 + 0.0*ksi2
    
            nu[9,1] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[9,2] = 0.0 + 0.0*ksi1 + 0.0*ksi2
            nu[9,3] = 0.0*ksi1 + 1.0*ksi2 + 0.0
    
            for s = 1 : DoF
                for t = 1 : 3
                  phi[s,t] = subs(nu[s,t], ksi1 => (x1+1)/2, ksi2 => (x2+1)/2)
                end
            end
            
            for i = 1 : DoF
                for j = 1 : Comp
                    Gradphi[i,j,1] = differentiate(phi[i,j],x1)
                    Gradphi[i,j,2] = differentiate(phi[i,j],x2)
                end
            end
    
            points[1,:] = [-1.0 -1.0]
            points[2,:] = [1.0 -1.0]
            points[3,:] = [-1.0 1.0]
    
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

        elseif Type == Grids.Quad() #Quad
            Comp = 3
            kp1 = k + 1
            DoF = (kp1 * kp1) * Comp
            Glob = KernelAbstractions.zeros(backend,Int,0,0)
            points = KernelAbstractions.zeros(backend,Float64,4,2)
            phi = Array{Polynomial,2}(undef,DoF,Comp)
            Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
            L1 = Array{Polynomial,1}(undef,kp1)
            L2 = Array{Polynomial,1}(undef,kp1)

            @polyvar x[1:2] 
            xw,_= gausslobatto(kp1)
            for i = 1 : kp1
                L1[i] = Lagrange(x[1],xw,i)   
                L2[i] = Lagrange(x[2],xw,i)  
            end  
            iDoF = 1
            for j = 1 : kp1
                for i = 1 : kp1
                    phi[iDoF,1] = L1[i] * L2[j]
                    phi[iDoF,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
                    phi[iDoF,3] = 0.0 + 0.0*x[1] + 0.0*x[2]
                    iDoF += 1
               
                    phi[iDoF,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
                    phi[iDoF,2] = L1[i] * L2[j]
                    phi[iDoF,3] = 0.0 + 0.0*x[1] + 0.0*x[2]
                    iDoF += 1
                
                    phi[iDoF,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
                    phi[iDoF,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
                    phi[iDoF,3] = L1[i] * L2[j]
                    iDoF += 1
                end  
            end 
            
            iDoF = 1
            for j = 1 : kp1
                for i = 1 : kp1
                    points[iDoF,1] = xw[i]
                    points[iDoF,2] = xw[j]
                    iDoF += 1
                    
                end  
            end   
            for i = 1 : DoF
                for j = 1 : Comp
                    Gradphi[i,j,1] = differentiate(phi[i,j],x[1])
                    Gradphi[i,j,2] = differentiate(phi[i,j],x[2])
                end
            end

            Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
            GlobCPU = zeros(Int,DoF,Grid.NumFaces)
            NumG = DoF * Grid.NumFaces
            NumI = DoF * Grid.NumFaces
            for iF = 1 : Grid.NumFaces
                for iD = 1 : DoF
                    GlobCPU[iD,iF] = DoF * (Grid.Faces[iF].F - 1 ) + iD
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