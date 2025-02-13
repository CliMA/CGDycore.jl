#constructing the ND elements for the Triangular grid
function ConstructND(k,ElemType::Grids.Tri)
  s = @polyvar x[1:2]

  P_k = Polynomial_k(k,s)
  lP_k = length(P_k)
  if k > 0
    P_km1 = Polynomial_k(k-1,s)
    lP_km1 = length(P_km1)
  else
    lP_km1 = 0  
  end  
  H_km1 = HomegenuousPolynomial(k,s)
  lH_km1 = length(H_km1)
  DoF = 2 * lP_k + lH_km1
  DoFE = lH_km1
  DoFF = DoF - 3 * DoFE
  phi = Array{Polynomial,2}(undef,DoF,2)
  phiB = Array{Polynomial,2}(undef,DoF,2)
  rounded_poly = Array{Polynomial,2}(undef,DoF,2)
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  rounded_Curlphi = Array{Polynomial,2}(undef,DoF,1) 
  iDoF = 1 
  @inbounds for i = 1 : lP_k
    phi[iDoF,1] = P_k[i]  
    phi[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
    iDoF += 1
    phi[iDoF,2] = P_k[i]  
    phi[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
    iDoF += 1
  end  
  @inbounds for i = 1 : lH_km1
    phi[iDoF,1] = -H_km1[i] * x[2]
    phi[iDoF,2] = H_km1[i] * x[1]
    iDoF += 1
  end  
  @polyvar t
  phiL = CGLine(k,t)
  QuadOrd = 3
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  I = zeros(DoF,DoF)
  rDoF = 1
# Compute functional over edges
  # Edge 1 (-1,-1) -> (1,-1)
  @inbounds for iDoF = 1 : DoF
    phiE1 = subs(phi[iDoF,1], x[1] => t, x[2] => -1.0)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE1(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
  end
  rDoF += k + 1
  # Edge 2 (1,-1) -> (-1,1)
  @inbounds for iDoF = 1 : DoF
    phiE1 = subs(phi[iDoF,1], x[1] => -t, x[2] => t)
    phiE2 = subs(phi[iDoF,2], x[1] => -t, x[2] => t)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += -0.5 * (phiE1(PointsL[iQ]) - phiE2(PointsL[iQ])) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ] 
      end
    end
  end
  rDoF += k + 1
# Edge 3 (-1,1) -> (-1,-1)
  @inbounds for iDoF = 1 : DoF
    phiE2 = subs(phi[iDoF,2], x[1] => -1, x[2] => -t)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE2(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]  
      end  
    end  
  end  
  rDoF += k + 1
  NumQuadT, WeightsT, PointsT = FEMSei.QuadRule(Grids.Tri(),QuadOrd)
# Interior  
  @inbounds for i = 1 : lP_km1
    @inbounds for iDoF = 1 : DoF
      @inbounds for iQ = 1 : NumQuadT
        Fac = P_km1[i](PointsT[iQ,1],PointsT[iQ,2])  
        I[rDoF,iDoF] += 0.25 * Fac * phi[iDoF,1](PointsT[iQ,1],PointsT[iQ,2]) * WeightsT[iQ] 
        I[rDoF+1,iDoF] += 0.25 * Fac * phi[iDoF,2](PointsT[iQ,1],PointsT[iQ,2]) * WeightsT[iQ]
      end
    end
    rDoF += 2
  end

  @inbounds for iDoF = 1 : DoF  
    @inbounds for jDoF = 1 : DoF  
      if abs(I[iDoF,jDoF]) < 1.e-12
        I[iDoF,jDoF] = 0
      end
    end
  end  
  r = zeros(DoF)
  @inbounds for iDoF = 1 : DoF  
    r[iDoF] = 1
    c = I \ r
    phiB[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
    phiB[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
    @inbounds for jDoF = 1 : DoF  
      phiB[iDoF,:] += c[jDoF] * phi[jDoF,:]
    end  
    phiB[iDoF,1] = round.(phiB[iDoF,1], digits=5)
    phiB[iDoF,2] = round.(phiB[iDoF,2], digits=5)
    r[iDoF] = 0
  end  
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  @inbounds for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phiB[i,1],x[2]) + differentiate(phiB[i,2],x[1])
  end
  return DoF, DoFE, DoFF, phiB, Curlphi
end

#constructing the ND elements for the Quadrilateral grid
function ConstructND(k,ElemType::Grids.Quad)

  s = @polyvar x[1:2]
  P_kp1x1 = Polynomial_1D(k+1,s,1)
  P_kx1 = Polynomial_1D(k,s,1)
  P_kp1x2 = Polynomial_1D(k+1,s,2)
  P_kx2 = Polynomial_1D(k,s,2)
  if k > 0
    P_km1x2 = Polynomial_1D(k-1,s,2)
    P_km1x1 = Polynomial_1D(k-1,s,1)
  end  
  
  DoF = 2 * (k+2) * (k+1)
  DoFE = k + 1
  DoFF = DoF - 4 * DoFE

  phi = Array{Polynomial,2}(undef,DoF,2)
  phiB = Array{Polynomial,2}(undef,DoF,2)
  rounded_poly = Array{Polynomial,2}(undef,DoF,2)
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  rounded_Curlphi = Array{Polynomial,2}(undef,DoF,1)
  iDoF = 1 
  @inbounds for i = 1 : k+2
    @inbounds for j = 1 : k+1
      phi[iDoF,2] = P_kp1x1[i] * P_kx2[j] 
      phi[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
      iDoF += 1
      phi[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
      phi[iDoF,1] = P_kp1x2[i] * P_kx1[j] 
      iDoF += 1
    end
  end
  @polyvar t
  phiL = CGLine(k,t)
  QuadOrd = 3
  NumQuadL, WeightsL, PointsL = FEMSei.QuadRule(Grids.Line(),QuadOrd)
  I = zeros(DoF,DoF)
  rDoF = 1
# Compute functional over edges
  # Edge 1 (-1,-1) -> (1,-1)
  @inbounds for iDoF = 1 : DoF
    phiE1 = subs(phi[iDoF,1], x[1] => t, x[2] => -1.0)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE1(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
  end
  rDoF += k + 1
  # Edge 2 (1,-1) -> (1,1)
  @inbounds for iDoF = 1 : DoF
    phiE2 = subs(phi[iDoF,2], x[1] => 1.0, x[2] => t)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE2(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]  
      end  
    end  
  end 
  rDoF += k + 1
  # Edge 3 (-1,1) -> (1,1)
  @inbounds for iDoF = 1 : DoF
    phiE1 = subs(phi[iDoF,1], x[1] => t, x[2] => 1.0)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE1(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]
      end
    end
  end
  rDoF += k + 1
  # Edge 4 (-1,-1) -> (-1,1)
  @inbounds for iDoF = 1 : DoF
    phiE2 = subs(phi[iDoF,2], x[1] => -1.0, x[2] => t)
    @inbounds for i = 0 : k
      @inbounds for iQ = 1 : NumQuadL
        I[rDoF+i,iDoF] += 0.5 * phiE2(PointsL[iQ]) * phiL[i+1](PointsL[iQ]) * WeightsL[iQ]  
      end  
    end  
  end  
  rDoF += k + 1
  NumQuadT, WeightsT, PointsT = FEMSei.QuadRule(Grids.Quad(),QuadOrd)
# Interior  
  @inbounds for i = 1 : k+1
    @inbounds for j = 1 : k
      @inbounds for iDoF = 1 : DoF
        phiI1 = phi[iDoF,1]  
        phiI2 = phi[iDoF,2]  
        @inbounds for iQ = 1 : NumQuadT
          I[rDoF,iDoF] += 0.25 * phiI2(PointsT[iQ,1],PointsT[iQ,2]) * 
          P_km1x1[j](PointsT[iQ,1],PointsT[iQ,2]) * P_kx2[i](PointsT[iQ,1],PointsT[iQ,2]) * 
          WeightsT[iQ]
          I[rDoF+1,iDoF] += 0.25 * phiI1(PointsT[iQ,1],PointsT[iQ,2]) * 
          P_kx1[i](PointsT[iQ,1],PointsT[iQ,2]) * P_km1x2[j](PointsT[iQ,1],PointsT[iQ,2]) * 
          WeightsT[iQ] 
        end
      end
      rDoF += 2
    end
  end
  @inbounds for iDoF = 1 : DoF  
    @inbounds for jDoF = 1 : DoF  
      if abs(I[iDoF,jDoF]) < 1.e-12
        I[iDoF,jDoF] = 0
      end
    end
  end  
  r = zeros(DoF)
  @inbounds for iDoF = 1 : DoF  
    r[iDoF] = 1
    c = I \ r
    phiB[iDoF,1] = 0.0 * x[1] + 0.0 * x[2]
    phiB[iDoF,2] = 0.0 * x[1] + 0.0 * x[2]
    @inbounds for jDoF = 1 : DoF  
      phiB[iDoF,:] += c[jDoF] * phi[jDoF,:]
    end  
    phiB[iDoF,1] = round.(phiB[iDoF,1], digits=5)
    phiB[iDoF,2] = round.(phiB[iDoF,2], digits=5)
    r[iDoF] = 0
  end  
  Curlphi = Array{Polynomial,2}(undef,DoF,1)
  @inbounds for i = 1 : DoF
    Curlphi[i,1] = -differentiate(phiB[i,1],x[2]) + differentiate(phiB[i,2],x[1])
  end
  return DoF, DoFE, DoFF, phiB, Curlphi
end

mutable struct NDStruct{FT<:AbstractFloat,
                        IT2<:AbstractArray} <: HCurlConfElement
  Order::Int                    
  Glob::IT2
  DoF::Int
  Comp::Int                      
  phi::Array{Polynomial,2}  
  Curlphi::Array{Polynomial,2}                       
  NumG::Int
  NumI::Int
  Type::Grids.ElementType
  M::AbstractSparseMatrix
  LUM::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64}
end

function NDStruct{FT}(backend,k,ElemType::Grids.ElementType,Grid) where FT<:AbstractFloat
  @polyvar x[1:2]
  Glob = KernelAbstractions.zeros(backend,Int,0,0)
  DoF, DoFE, DoFF, phi, Curlphi = FEMSei.ConstructND(k,ElemType)
  Comp = 2
  Glob = KernelAbstractions.zeros(backend,Int,DoF,Grid.NumFaces)
  GlobCPU = zeros(Int,DoF,Grid.NumFaces)
  NumG = Grid.NumEdges * DoFE + Grid.NumFaces * DoFF
  NumI = NumG
  if ElemType == Grids.Tri
    @inbounds for iF = 1 : Grid.NumFaces
      iGlob = 1
      @inbounds for i = 1 : length(Grid.Faces[iF].E)
        iE = Grid.Faces[iF].E[i]
        OrientE = Grid.Faces[iF].OrientE[i]
        @inbounds for j = 1 : DoFE
          if OrientE > 0
            GlobCPU[iGlob,iF] = DoFE * (Grid.Edges[iE].E - 1) + j
          else
            GlobCPU[iGlob,iF] = DoFE * (Grid.Edges[iE].E - 1) + DoFE - j + 1
          end
          iGlob += 1
        end
      end
      @inbounds for j = 1 : DoFF
        GlobCPU[iGlob,iF] = DoFE * Grid.NumEdges + DoFF * (Grid.Faces[iF].F - 1) + j
        iGlob += 1
      end
    end
  else
    @inbounds for iF = 1 : Grid.NumFaces
      iGlob = 1
      @inbounds for i = 1 : length(Grid.Faces[iF].E)
        iE = Grid.Faces[iF].E[i]
        @inbounds for j = 1 : DoFE
          GlobCPU[iGlob,iF] = DoFE * (Grid.Edges[iE].E - 1) + j
          iGlob += 1
        end
      end
      @inbounds for j = 1 : DoFF
        GlobCPU[iGlob,iF] = DoFE * Grid.NumEdges + DoFF * (Grid.Faces[iF].F - 1) + j
        iGlob += 1
      end
    end
  end
  copyto!(Glob,GlobCPU)
  M = sparse([1],[1],[1.0])
  LUM = lu(M)
  Order = k
  return NDStruct{FT,
                  typeof(Glob)}( 
    Order,              
    Glob,
    DoF,
    Comp,
    phi,
    Curlphi,                      
    NumG,
    NumI,
    ElemType,
    M,
    LUM,
      )
end
