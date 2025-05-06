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
function VecDGStruct{FT}(backend,k::Int,Type::Grids.Quad,Grid) where FT<:AbstractFloat
    
  if k == 0 #order 0
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
    @inbounds for i = 1 : DoF
      @inbounds for j = 1 : Comp
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
    @inbounds for iF = 1 : Grid.NumFaces
      @inbounds for iDoF = 1 : DoF
        GlobCPU[iDoF,iF] = (Grid.Faces[iF].F - 1 ) * DoF + iDoF
      end
    end
  else 
    Comp = 3
    kp1 = k + 1
    DoF = (kp1 * kp1) * Comp
    Glob = KernelAbstractions.zeros(backend,Int,0,0)
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    L1 = Array{Polynomial,1}(undef,kp1)
    L2 = Array{Polynomial,1}(undef,kp1)

    @polyvar x[1:2] 
    xw,_= gausslobatto(kp1)
    @inbounds for i = 1 : kp1
      L1[i] = Lagrange(x[1],xw,i)   
      L2[i] = Lagrange(x[2],xw,i)  
    end  
    iDoF = 1
    @inbounds for j = 1 : kp1
      @inbounds for i = 1 : kp1
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
           
    points = KernelAbstractions.zeros(backend,Float64,kp1*kp1,2)
    iDoF = 1
    @inbounds for j = 1 : kp1
      @inbounds for i = 1 : kp1
        points[iDoF,1] = xw[i]
        points[iDoF,2] = xw[j]
        iDoF += 1
      end  
    end   
    @inbounds for iDoF = 1 : DoF
      @inbounds for iComp = 1 : Comp
        Gradphi[iDoF,iComp,1] = differentiate(phi[iDoF,iComp],x[1])
        Gradphi[iDoF,iComp,2] = differentiate(phi[iDoF,iComp],x[2])
      end
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = DoF * Grid.NumFaces
  NumI = DoF * Grid.NumFaces
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iD = 1 : DoF
      GlobCPU[iD,iF] = DoF * (Grid.Faces[iF].F - 1 ) + iD
    end
  end
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M) 

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

function VecDGStruct{FT}(backend,k::Int,Type::Grids.Tri,Grid) where FT<:AbstractFloat
    
  @polyvar x[1:2]  
  if k == 0 #order 0
    DoF = 3 
    Comp = 3
    @polyvar x1 x2
    phi = Array{Polynomial,2}(undef,DoF,Comp)
    phi[1,1] = 1.0 + 0.0*x[1] + 0.0*x[2]
    phi[1,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
    phi[1,3] = 0.0 + 0.0*x[1] + 0.0*x[2]
    phi[2,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
    phi[2,2] = 1.0 + 0.0*x[1] + 0.0*x[2]
    phi[2,3] = 0.0 + 0.0*x[1] + 0.0*x[2]
    phi[3,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
    phi[3,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
    phi[3,3] = 1.0 + 0.0*x[1] + 0.0*x[2]
    Gradphi = Array{Polynomial,3}(undef,DoF,Comp,2)
    @inbounds for i = 1 : DoF
      @inbounds for j = 1 : Comp
        Gradphi[i,j,1] = differentiate(phi[i,j],x[1])
        Gradphi[i,j,2] = differentiate(phi[i,j],x[2])
      end
    end
        
    points = KernelAbstractions.zeros(backend,Float64,DoF,2)
    points[1,:] = [0.0 0.0]
          
  else
    Comp = 3 
    DoF, DoFE, DoFF, phiDG, GradphiDG, points = FEMSei.ConstructCG(k,Type)
    phi = Array{Polynomial,2}(undef,Comp*DoF,Comp)
    Gradphi = Array{Polynomial,3}(undef,Comp*DoF,Comp,2)

    for iDoF = 1 : DoF
      phi[Comp*iDoF-2,1] = phiDG[iDoF]
      Gradphi[Comp*iDoF-2,1,1] = GradphiDG[iDoF,1,1]
      Gradphi[Comp*iDoF-2,1,2] = GradphiDG[iDoF,1,2]
      phi[Comp*iDoF-2,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-2,2,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-2,2,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
      phi[Comp*iDoF-2,3] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-2,3,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-2,3,2] = 0.0 + 0.0*x[1] + 0.0*x[2]

      phi[Comp*iDoF-1,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-1,1,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-1,1,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
      phi[Comp*iDoF-1,2] = phiDG[iDoF]
      Gradphi[Comp*iDoF-1,2,1] = GradphiDG[iDoF,1,1]
      Gradphi[Comp*iDoF-1,2,2] = GradphiDG[iDoF,1,2]
      phi[Comp*iDoF-1,3] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-1,3,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF-1,3,2] = 0.0 + 0.0*x[1] + 0.0*x[2]

      phi[Comp*iDoF,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF,1,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF,1,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
      phi[Comp*iDoF,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF,2,1] = 0.0 + 0.0*x[1] + 0.0*x[2]
      Gradphi[Comp*iDoF,2,2] = 0.0 + 0.0*x[1] + 0.0*x[2]
      phi[Comp*iDoF,3] = phiDG[iDoF]
      Gradphi[Comp*iDoF,3,1] = GradphiDG[iDoF,1,1]
      Gradphi[Comp*iDoF,3,2] = GradphiDG[iDoF,1,2]
    end  
    DoF *= Comp
  end  

  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = DoF * Grid.NumFaces
  NumI = DoF * Grid.NumFaces
  @inbounds for iF = 1 : Grid.NumFaces
    @inbounds for iDoF = 1 : DoF
      GlobCPU[iDoF,iF] = (Grid.Faces[iF].F - 1 ) * DoF + iDoF
    end
  end

  Glob = KernelAbstractions.zeros(backend,Int,size(GlobCPU))
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)

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
