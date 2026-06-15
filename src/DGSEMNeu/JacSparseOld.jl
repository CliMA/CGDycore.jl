import CGDycore:
  Parameters, Thermodynamics, Examples, Sources, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGSEM, CGSEM, DyCore, IMEXRosenbrock

using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using KernelAbstractions
using KernelAbstractions: @atomic, @atomicswap, @atomicreplace
using KernelAbstractions.Extras: @unroll
using MPI

abstract type RiemannSolver end

Base.@kwdef struct RiemannLMARSLinFast <: RiemannSolver end

function (::RiemannLMARSLinFast)(Phys,RhoPos,uPos,vPos,wPos,RhoThPos,dpdRhoThPos,ThPos)
  @inline function RiemannByLMARSSemi!(F,VLL,VRR,AuxL,AuxR,n1,n2,n3)

    FT = eltype(F)
    pLL = AuxL[dpdRhoThPos] * VLL[RhoThPos]
    pRR = AuxR[dpdRhoThPos] * VRR[RhoThPos]
    vLL = (VLL[uPos] * n1 + VLL[vPos] * n2 + VLL[wPos] * n3)
    vRR = (VRR[uPos] * n1 + VRR[vPos] * n2 + VRR[wPos] * n3)
    pM = FT(0.5) * ((pLL + pRR) - Phys.cS * (vRR - vLL))
    vM = FT(0.5) * ((vRR + vLL) - Phys.invcS * (pRR - pLL))
    F[uPos] = n1 * pM
    F[vPos] = n2 * pM
    F[wPos] = n3 * pM
    F[RhoPos] = vM
    F[RhoThPos] = FT(0.5) * F[RhoPos] * (AuxL[ThPos] + AuxR[ThPos])
  end
  return RiemannByLMARSSemi!
end


abstract type AverageFlux end

Base.@kwdef struct KennedyGruberGravLinFast <: AverageFlux end

function (::KennedyGruberGravLinFast)(RhoPos,uPos,vPos,wPos,ThPos,pRhoThPos,ThAuxPos,GPPos)
  @inline function FluxNonLinAverSemi!(flux,
      VLoc, AuxLoc, dXdxILoc,
      K1, I1, J1,   # left state indices  (localidx into @localmem)
      K2, I2, J2,   # right state indices
      ::Val{dir}) where {dir}

    FT = eltype(flux)
    m_L1  = dXdxILoc[dir, 1, K1, I1, J1]
    m_L2  = dXdxILoc[dir, 2, K1, I1, J1]
    m_L3  = dXdxILoc[dir, 3, K1, I1, J1]
    m_R1  = dXdxILoc[dir, 1, K2, I2, J2]
    m_R2  = dXdxILoc[dir, 2, K2, I2, J2]
    m_R3  = dXdxILoc[dir, 3, K2, I2, J2]
    diffGP = AuxLoc[K2, I2, J2, GPPos] - AuxLoc[K1, I1, J1, GPPos]
    pAvL = AuxLoc[K1, I1, J1, pRhoThPos] * VLoc[K1, I1, J1, ThPos]  + 
      FT(0.5) * VLoc[K1, I1, J1, RhoPos] * diffGP
    pAvR = AuxLoc[K2, I2, J2, pRhoThPos] * VLoc[K2, I2, J2, ThPos]  + 
      FT(0.5) * VLoc[K2, I2, J2, RhoPos] * diffGP
    qHatL = m_L1 * VLoc[K1, I1, J1, uPos] + 
      m_L2 * VLoc[K1, I1, J1, vPos] + m_L3 * VLoc[K1, I1, J1, wPos]
    qHatR = m_R1 * VLoc[K2, I2, J2, uPos] + 
      m_R2 * VLoc[K2, I2, J2, vPos] + m_R3 * VLoc[K2, I2, J2, wPos]
    flux[1] = FT(0.5) * (qHatL + qHatR)
    flux[2] = FT(0.5) * (m_L1 * pAvL + m_R1 * pAvR)
    flux[3] = FT(0.5) * (m_L2 * pAvL + m_R2 * pAvR)
    flux[4] = FT(0.5) * (m_L3 * pAvL + m_R3 * pAvR)
    flux[5] = FT(0.5) * (qHatL * AuxLoc[K1, I1, J1, ThAuxPos] + qHatR * AuxLoc[K2, I2, J2, ThAuxPos])

  end
  return FluxNonLinAverSemi!
end

@inline function Index(k,i,j,iz,iF,Pos,M,nz,N,NF)
  ind = k + (iz - 1) * M + (i - 1) * M * nz + (j-1) * M * nz * N + 
    (iF - 1) * M * nz * N * N + (Pos - 1) * M * nz * N * N * NF
end  

@inline function Index(k,ijF,iz,Pos,M,nz,N,NF)
  ind = k + (iz - 1) * M + (ijF - 1) * M * nz + 
    (Pos - 1) * M * nz * N * N * NF
end    

@inline function Insert!(iRow,iCol,v,RowInd,ColInd,Val)
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,v)
end  

function JacFluxVolumeSparse(DG,Metric,Aux,Phys)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  RowIndVR = Int[]
  ColIndVR = Int[]
  ValVR = Float64[]
  RowIndHR = Int[]
  ColIndHR = Int[]
  ValHR = Float64[]

  DVT = deepcopy(DG.DVT)
  N = size(DVT,1)
  DVZT = deepcopy(DG.DVZT)
  M = size(DVZT,1)
  DVT = deepcopy(DG.DVT)
  for i = 1 : N
    DVT[i,i] = 0.0  
    DVT[i,i] = sum(DVT[:,i])
  end  
  DVZT = deepcopy(DG.DVZT)
  for i = 1 : M
    DVZT[i,i] = 0.0  
    DVZT[i,i] = sum(DVZT[:,i])
  end  
  dXdxIG = Metric.dXdxI
  NH = Metric.NH
  VolSurfH = Metric.VolSurfH
  NVG = Metric.NV
  VolSurfVG = Metric.VolSurfV
  nz = size(dXdxIG,5)
  NF = size(dXdxIG,6)
  NE = size(DG.GlobE,2)
  dXdxI = reshape(dXdxIG,3,3,M,N,N,nz,NF)
  VolSurfV = reshape(VolSurfVG,nz+1,N,N,NF)
  NV = reshape(NVG,nz+1,N,N,NF,3)
  Theta = reshape(Aux[:,:,:,4],M,nz,N,N,NF)
  pRhoTh = reshape(Aux[:,:,:,3],M,nz,N,N,NF)
  Mnz = nz * M
  MnzNN = N * N * M * nz
  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N  
        for iz = 1 : nz  
          for k = 1 : M  
            # Rho RhoTheta 
            iRowRho = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            iRowRhoTh = Index(k,i,j,iz,iF,RhoThPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,uPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,1,k,l,j,iz,iF] * DVT[l,i]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,l,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,l,j,iz,iF,vPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,2,k,l,j,iz,iF] * DVT[l,i]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,l,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,l,j,iz,iF,wPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,3,k,l,j,iz,iF] * DVT[l,i]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,l,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
            end  
#           j direction 
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,uPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,1,k,i,l,iz,iF] * DVT[l,j]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,i,l,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,i,l,iz,iF,vPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,2,k,i,l,iz,iF] * DVT[l,j]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,i,l,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,i,l,iz,iF,wPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,3,k,i,l,iz,iF] * DVT[l,j]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,i,l,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
            end      
#           k direction 
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,uPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,1,l,i,j,iz,iF] * DVZT[l,k]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[l,iz,i,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(l,i,j,iz,iF,vPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,2,l,i,j,iz,iF] * DVZT[l,k]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[l,iz,i,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(l,i,j,iz,iF,wPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,3,l,i,j,iz,iF] * DVZT[l,k]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[l,iz,i,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
            end         
            # u  
            iRow = Index(k,i,j,iz,iF,uPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,1,k,l,j,iz,iF] * pRhoTh[k,iz,l,j,iF] * DVT[l,i]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end  
#           j direction
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,1,k,i,l,iz,iF] * pRhoTh[k,iz,i,l,iF] * DVT[l,j]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end  
#           k direction
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,1,l,i,j,iz,iF] * pRhoTh[l,iz,i,j,iF] * DVZT[l,k]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end  

            # v  
            iRow = Index(k,i,j,iz,iF,vPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,2,k,l,j,iz,iF] * pRhoTh[k,iz,l,j,iF] * DVT[l,i]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end 
#           j direction
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,2,k,i,l,iz,iF] * pRhoTh[k,iz,i,l,iF] * DVT[l,j]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end 
#           k direction
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,2,l,i,j,iz,iF] * pRhoTh[l,iz,i,j,iF] * DVZT[l,k]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end         

            # w  
            iRow = Index(k,i,j,iz,iF,wPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,3,k,l,j,iz,iF] * pRhoTh[k,iz,l,j,iF] * DVT[l,i]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end
#           j direction
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,3,k,i,l,iz,iF] * pRhoTh[k,iz,i,l,iF] * DVT[l,j]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end
#           k direction
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,3,l,i,j,iz,iF] * pRhoTh[l,iz,i,j,iF] * DVZT[l,k]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end        
          end        
        end        
      end        
    end        
  end        
# Adding the surface terms
# Vertical direction  
  @show NF,N,nz
  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N
        for iz = 1 : nz + 1
          n1 = NV[iz,i,j,iF,1]
          n2 = NV[iz,i,j,iF,2]
          n3 = NV[iz,i,j,iF,3]  
          Surf = VolSurfV[iz,i,j,iF]
          if iz == 1
            # u
            iRowT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            # u
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n1 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n2 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n3 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n1 * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)

            # v
            iRowT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            # u
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n1 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n2 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n3 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n2 * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)

            # w
            iRowT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            # u
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n1 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n2 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n3 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n3 * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
          elseif iz == nz + 1
            # u
            iRowB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n1 * pRhoTh[M,iz-1,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)

            # v
            iRowB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n2 * pRhoTh[M,iz-1,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)

            # w
            iRowB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n3 * pRhoTh[M,iz-1,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
          else
            # Rho 
            iRowB = Index(M,i,j,iz-1,iF,RhoPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,RhoPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * Phys.invcS * pRhoTh[M,iz-1,i,j,iF] * Surf
            ValLocT = 0.5 * Phys.invcS * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLocT,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLocT,RowIndVR,ColIndVR,ValVR)


            # RhoTh 
            iRowB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.25 * n1 * Surf * (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF])
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.25 * n2 * Surf * (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF])
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.25 * n3 * Surf * (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF])
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.25 * Phys.invcS * pRhoTh[M,iz-1,i,j,iF] * 
              (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF]) * Surf
            ValLocT = 0.25 * Phys.invcS * pRhoTh[1,iz,i,j,iF] * 
              (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF]) * Surf
            Insert!(iRowB,iColB,-ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLocT,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLocT,RowIndVR,ColIndVR,ValVR)

            # u 
            iRowB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * n1 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * n2 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * n3 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * n1 * pRhoTh[M,iz-1,i,j,iF] * Surf 
            ValLocT = 0.5 * n1 * pRhoTh[1,iz,i,j,iF] * Surf 
            Insert!(iRowB,iColB,-ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLocT,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLocT,RowIndVR,ColIndVR,ValVR)

            # v 
            iRowB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * n1 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * n2 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * n3 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * n2 * pRhoTh[M,iz-1,i,j,iF] * Surf 
            ValLocT = 0.5 * n2 * pRhoTh[1,iz,i,j,iF] * Surf 
            Insert!(iRowB,iColB,-ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLocT,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLocT,RowIndVR,ColIndVR,ValVR)

            # w 
            iRowB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * n1 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * n2 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * n3 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLoc,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,-ValLoc,RowIndVR,ColIndVR,ValVR)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * n3 * pRhoTh[M,iz-1,i,j,iF] * Surf 
            ValLocT = 0.5 * n3 * pRhoTh[1,iz,i,j,iF] * Surf 
            Insert!(iRowB,iColB,-ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowB,iColT,-ValLocT,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColB,ValLocB,RowIndVR,ColIndVR,ValVR)
            Insert!(iRowT,iColT,ValLocT,RowIndVR,ColIndVR,ValVR)
          end
        end
      end
    end
  end  
  @views Theta = Aux[:,:,:,4]
  @views pRhoTh = Aux[:,:,:,3]
  for iE = 1 : NE
    iFL = Grid.EF[1,iE]
    iFR = Grid.EF[2,iE]
    for i = 1 : N
      indL = DG.GlobE[1,i,iE]
      indR = DG.GlobE[2,i,iE]
      for iz = 1 : nz
        for k = 1 : M  
          n1 = NH[k,iz,i,iE,1]
          n2 = NH[k,iz,i,iE,2]
          n3 = NH[k,iz,i,iE,3]  
          Surf = VolSurfH[k,iz,i,iE]

          # Rho
          iRowL = Index(k,indL,iz,RhoPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,RhoPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * Surf
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          # v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * Surf
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          # w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * Surf
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * Phys.invcS * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * Phys.invcS * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLocR,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLocL,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLocR,RowIndHR,ColIndHR,ValHR)

          # RhoTh
          iRowL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.25 * n1 * Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          # v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.25 * n2 * Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          # w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.25 * n3 * Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.25 * Phys.invcS * pRhoTh[k,iz,indL] * 
            Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          ValLocR = 0.25 * Phys.invcS * pRhoTh[k,iz,indR] * 
            Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLocL,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLocR,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLocL,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLocR,RowIndHR,ColIndHR,ValHR)

          # u
          iRowL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,uPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * n1 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          #v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * n2 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          #w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * n3 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * n1 * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * n1 * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowL,iColR,-ValLocR,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowR,iColL,ValLocL,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowR,iColR,ValLocR,RowIndVR,ColIndVR,ValVR)

          # v
          iRowL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,vPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * n1 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          #v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * n2 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          #w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * n3 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * n2 * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * n2 * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowL,iColR,-ValLocR,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowR,iColL,ValLocL,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowR,iColR,ValLocR,RowIndVR,ColIndVR,ValVR)

          # w
          iRowL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,wPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * n1 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          #v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * n2 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          #w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * n3 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowL,iColR,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColL,ValLoc,RowIndHR,ColIndHR,ValHR)
          Insert!(iRowR,iColR,-ValLoc,RowIndHR,ColIndHR,ValHR)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * n3 * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * n3 * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowL,iColR,-ValLocR,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowR,iColL,ValLocL,RowIndVR,ColIndVR,ValVR)
          Insert!(iRowR,iColR,ValLocR,RowIndVR,ColIndVR,ValVR)
        end
      end  
    end
  end  
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  JacVR = sparse(RowIndVR, ColIndVR, ValVR, NZe, NZe)
  JacHR = sparse(RowIndHR, ColIndHR, ValHR, NZe, NZe)
  return Jac,JacHR,JacVR
end

@kernel inbounds = true function RiemannNonLinV3Kernel!(RiemannSolver!,F,@Const(U),
  @Const(Aux),@Const(NV),@Const(VolSurfV),
  ::Val{NUMV}, ::Val{NAUX}) where {NUMV, NAUX}

  Iz,ind = @index(Global, NTuple)


  Nz = @uniform @ndrange()[1]
  NumI = @uniform @ndrange()[2]

  VLL = @private eltype(F) (NUMV,)
  VRR = @private eltype(F) (NUMV,)
  AuxL = @private eltype(F) (NAUX,)
  AuxR = @private eltype(F) (NAUX,)
  FLoc = @private eltype(F) (NUMV,)

  if ind <= NumI
    n1 = NV[Iz,ind,1]
    n2 = NV[Iz,ind,2]
    n3 = NV[Iz,ind,3]
    if Iz == 1
      @unroll for iAux = 1 : NAUX
        AuxL[iAux] = Aux[1,Iz,ind,iAux]
        AuxR[iAux] = AuxL[iAux]
      end
      @unroll for iv = 1 : NUMV
        VLL[iv] = U[1,Iz,ind,iv]
        VRR[iv] = VLL[iv]
      end
      t = eltype(F)(2) * (n1 * VLL[2] +
        n2 * VLL[3] +
        n3 * VLL[4])
      VLL[2] -= n1 * t
      VLL[3] -= n2 * t
      VLL[4] -= n3 * t
      RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]  
      @unroll for iv = 1 : NUMV  
        F[1,Iz,ind,iv] += FLoc[iv] * Surf
      end  
    elseif Iz == Nz 
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[end,Iz-1,ind,iAux]
        AuxR[iAux] = AuxL[iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[end,Iz-1,ind,iv]
        VRR[iv] = VLL[iv]
      end  
      t = eltype(F)(2) * (n1 * VRR[2] +
        n2 * VRR[3] +
        n3 * VRR[4])
      VRR[2] -= n1 * t
      VRR[3] -= n2 * t
      VRR[4] -= n3 * t
      RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]   
      @unroll for iv = 1 : NUMV  
        F[end,Iz-1,ind,iv] -= FLoc[iv] * Surf
      end  
    else
      @unroll for iAux = 1 : NAUX  
        AuxL[iAux] = Aux[end,Iz-1,ind,iAux]
        AuxR[iAux] = Aux[1,Iz,ind,iAux]
      end  
      @unroll for iv = 1 : NUMV  
        VLL[iv] = U[end,Iz-1,ind,iv]
        VRR[iv] = U[1,Iz,ind,iv]
      end  
      RiemannSolver!(FLoc,VLL,VRR,AuxL,AuxR,n1,n2,n3)
      Surf = VolSurfV[Iz,ind]  
      @unroll for iv = 1 : NUMV  
        FLoc[iv] *= Surf
        F[end,Iz-1,ind,iv] -= FLoc[iv]
        F[1,Iz,ind,iv] += FLoc[iv]
      end  
    end    
  end  
end

@kernel inbounds = true function FluxSplitVolumeNonLinHVQuadKernel!(FluxAver!,
    F, @Const(V), @Const(Aux), @Const(dXdxI), @Const(DVT),
    @Const(DVZT), @Const(Glob),
    ::Val{N}, ::Val{M}, ::Val{NV}, ::Val{NAUX}) where {N, M, NV, NAUX}

  K, I, J      = @index(Local, NTuple)
  _, _, _, IZ, IF  = @index(Global, NTuple)

  NZ       = @uniform @ndrange()[4]

  VLoc     = @localmem eltype(F)      (M, N, N, NV)
  AuxLoc   = @localmem eltype(F)      (M, N, N, NAUX)
  # 2 directions × 3 components
  dXdxILoc = @localmem eltype(dXdxI)  (3, 3, M, N, N, 1)

  fTilde = @private eltype(F) (NV,)
  gTilde = @private eltype(F) (NV,)
  hTilde = @private eltype(F) (NV,)
  FLoc = @private eltype(F) (NV,)

  # ---- load phase ----
  ID = I + (J - 1) * N
  ind = Glob[ID, IF]
  @unroll for iaux = 1:NAUX
    AuxLoc[K, I, J, iaux] = Aux[K, IZ, ind, iaux]
  end
  @unroll for iv = 1:NV
    VLoc[K, I, J, iv] = V[K, IZ, ind, iv]
    FLoc[iv] = 0
  end
  @unroll for j = 1:3
    @unroll for i = 1:3
      dXdxILoc[i, j, K, I, J] = dXdxI[i, j, K, ID, IZ, IF]
    end
  end

  @synchronize

  # ---- compute phase ----
  @unroll for l = 1:N
    # x-direction: left=(I,J,K), right=(l,J,K)
    FluxAver!(fTilde,
      VLoc, AuxLoc, dXdxILoc,
      K, I, J,
      K, l, J,
      Val(1))
    # y-direction: left=(I,J,K), right=(I,l,K)
    FluxAver!(gTilde,
      VLoc, AuxLoc, dXdxILoc,
      K, I, J,
      K, I, l,
      Val(2))
    @unroll for iv = 1:NV
      FLoc[iv] += -DVT[l, I] * fTilde[iv] - DVT[l, J] * gTilde[iv]
    end
  end  
  @unroll for l = 1:M
    # z-direction: left=(I,J,K), right=(I,J,l)
    FluxAver!(hTilde,
      VLoc, AuxLoc, dXdxILoc,
      K, I, J,
      l, I, J,
      Val(3))
    @unroll for iv = 1:NV
      FLoc[iv] += -DVZT[l, K] * hTilde[iv]
    end
  end  
  ID = I + (J - 1) * N
  ind = Glob[ID, IF]
  @unroll for iv = 1:NV
    F[K, IZ, ind, iv] += FLoc[iv]
  end
end
function DLagrange(x,xP,iP)
  Df=eltype(x)(0)
  f=eltype(x)(1)
  for i = 1 : size(xP,1)
    if i ≠ iP
      Df = Df * (x-xP[i]) / (xP[iP] - xP[i]) + f / (xP[iP]-xP[i])
      f = f * (x - xP[i]) / (xP[iP] - xP[i])
    end
  end
  return Df
end

###### Main




  backend = CPU()
  FTB = Float64
  MPI.Init()
  comm = MPI.COMM_WORLD
  Proc = MPI.Comm_rank(comm) + 1
  ProcNumber = MPI.Comm_size(comm)
  ParallelCom = DyCore.ParallelComStruct()
  ParallelCom.Proc = Proc
  ParallelCom.ProcNumber  = ProcNumber

  Phys = DyCore.PhysParameters{FTB}()
  Model = DyCore.ModelStruct{FTB}()
  OrdPoly = 4
  OrdPolyZ = 4
  N = OrdPoly + 1
  M = OrdPolyZ + 1
  nz = 3
  nPanel = 5
  RefineLevel = 3
  ns = 10
  nLon = 10
  nLat = 10
  LatB = 80
  RadEarth = 1.0
  GridType = "CubedSphere"
  Decomp="EqualArea"
  Discretization="DG"
  DGMethod = "Kubatko2LGL"
  AdaptGridType = "Sleve"
  H = 1000.0
  OrdPrint = 1
  OrdPrintZ = 1
  TopoS = ""
  TopoProfile = Examples.Flat()()

  Grid, CellToProc = Grids.InitGridSphere(backend,FTB,OrdPoly,nz,nPanel,RefineLevel,ns,nLon,nLat,LatB,
    GridType,Decomp,RadEarth,Model,ParallelCom;Discretization=Discretization,ChangeOrient=2)

  Grid.AdaptGrid = Grids.AdaptGrid(FTB,AdaptGridType,FTB(H))
  Topography = (TopoS=TopoS,H=H,Rad=RadEarth)

  (DG, Metric, Exchange, Global) = DyCore.InitSphereDG(backend,FTB,OrdPoly,OrdPolyZ,DGMethod,
    OrdPrint,OrdPrintZ,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5


  NUMV = 5
  nz = Grid.nz
  NF = size(DG.Glob,2)
  F = zeros(M,nz,DG.NumI,NUMV)
  FR = zeros(M,nz,DG.NumI,NUMV)
  U = rand(M,nz,DG.NumI,NUMV)
  NAUX = 4
  Aux = rand(M,nz,DG.NumI,NAUX)
  @. Aux[:,:,:,1:2] = 0
  @. Aux[:,:,:,3] += 400
  @. Aux[:,:,:,4] += 300

  @views Jac,JacHR,JacVR = JacFluxVolumeSparse(DG,Metric,Aux,Phys)

  U1 = reshape(U,M*nz*DG.NumI*NUMV)
  F1V = Jac * U1
  F1 = reshape(F1V,M,nz,N,N,NF,NUMV)
  FR1V = JacVR * U1
  FR1 = reshape(FR1V,M,nz,N,N,NF,NUMV)

  RiemannSolverFast = RiemannLMARSLinFast()(Phys,RhoPos,uPos,vPos,
    wPos,RhoThPos,3,4)

# FluxAverageFast = DGSEM.KennedyGruberGravLinFast()(RhoPos,uPos,vPos,wPos,
#   RhoThPos,3,4,2,Grid.Type)
  FluxAverageFast = KennedyGruberGravLinFast()(RhoPos,uPos,vPos,wPos,
    RhoThPos,3,4,2)

  group = (nz+1,1)
  ndrange = (nz+1,DG.NumI)
  KRiemannNonLinV3Kernel! = RiemannNonLinV3Kernel!(backend,group)
  KRiemannNonLinV3Kernel!(RiemannSolverFast,FR,U,Aux,Metric.NV,
    Metric.VolSurfV,Val(NUMV),Val(NAUX);ndrange=ndrange)
  FR2 = reshape(FR,M,nz,N,N,NF,NUMV)

  group = (M,N,N,1,1)
  ndrange = (M,N,N,nz,NF)
  KFluxSplitVolumeNonLinHKernel! = FluxSplitVolumeNonLinHVQuadKernel!(backend,group)
  @views KFluxSplitVolumeNonLinHKernel!(FluxAverageFast,F,U,Aux,Metric.dXdxI,DG.DVT,DG.DVZT,DG.Glob,
    Val(N), Val(M), Val(NUMV), Val(NAUX);ndrange=ndrange)
  F2 = reshape(F,M,nz,N,N,NF,NUMV)
  aa = 3;

