import CGDycore:
  Thermodynamics, Sources, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DyCore, DGSEMNeu
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using LinearAlgebra
using SparseArrays
using BandedMatrices
using FillArrays

mutable struct Interior 
  M::Int
  nz::Int
  fac::Float64
  FacGrav::Float64
  A13::Array{Float64, 4}
  A23::Array{Float64, 4}
  A32::Array{Float64, 4}
  B1m_34::Array{Float64, 3}
  B1_1::Array{Float64, 2}
  B1_23::Array{Float64, 4}
  B1_4::Array{Float64, 2}
  B2_23::Array{Float64, 4}
  B3_14::Array{Float64, 4}
  B1p_12::Array{Float64, 3}
  C23_2::Array{Float64, 4}
  C14_3::Array{Float64, 4}
  luSA::Array{LU, 2}
  SchurBand::Array{BandedMatrix, 1}
end  

function ldivBand!(ABanded,b)

# Ly = b
  A = ABanded.data
  n = size(A,2)
  l = ABanded.l
  u = ABanded.u
  up1 = u + 1 
  uplp1 = up1 + l

  for k = 1 : n - 1
    for i = k + 1 : min(k+l,n)
      b[i] -= A[i - k + up1,k] * b[k]  
    end    
  end
# Ux = y
  for k = n : -1 : 1
    b[k] /= A[up1,k]
    for i = max(k - u, 1) : k - 1
      b[i] -= A[up1 + i - k,k] * b[k]  
    end    
  end
end

function luBand!(ABanded)
#1  *   *   *  a14  ...    
#2  *   *  a13 a24  ... 
#3  *  a12 a23 a34  ... 
#4 a11 a22 a33 a44  ...  
#5 a21 a32 a43 a54  ...  
#6 a31 a42 a53 a64  ...  
#7 a41 a52 a63 a74  ...  

#1  an-6n-3 an-5n-2 an-4n-1 an-3n
#2  an-5n-3 an-4n-2 an-3n-1 an-2n
#3  an-4n-3 an-3n-2 an-2n-1 an-1n
#4  an-3n-3 an-2n-2 an-1n-1 an  n
#5  an-2n-3 an-1n-2 an  n-1   *
#6  an-1n-3 an  n-2   *       *
#7  an  n-3    *      *       *


# a11 a12 a13 a14 ....
# a21 a22 a23 a24 a25 ...
# a31 a32 a33 a34 a35 a36 ...
# a41 a42 a43 a44 a45 a46 a47 ...
# ... a52 a53 a54 a55 a56 a57 a58 ...

#  an-1n-3 an-1n-2 an-1n-1 an-1n 
#    ann-3   ann-2   ann-1   ann

# a21 => (5,1)
# a31 => (6,1)
# a41 => (7,1)

# a12 => (3,2)
# a13 => (2,3) 
# a14 => (1,4)
  A = ABanded.data
  n = size(A,2)
  u = ABanded.u
  l = ABanded.l
  up1 = u + 1
  up2 = u + 2
  uplp1 = u + l + 1
  @inbounds for k = 1 : n -1 
    @inbounds for i = up2 : uplp1 
      A[i,k] /= A[up1,k]  
    end  
    @inbounds for i = up2 : uplp1 
      @inbounds for j = k + 1 : min(k+u,n)
        A[i+k-j,j] -= A[i,k] * A[k + up1 - j,j]
      end
    end  
  end    
end

function Interior(M,nz,NumG)
  M2 = M - 2
  fac = 0
  FacGrav = 0
  A13 = zeros(M,M2,nz,NumG)
  A23 = zeros(M2,M2,nz,NumG)
  A32 = zeros(M2,M2,nz,NumG)
  B1m_34 = zeros(2,nz,NumG)
  B1_1 = zeros(nz,NumG)
  B1_23 = zeros(M,2,nz,NumG)
  B1_4 = zeros(nz,NumG)
  B2_23 = zeros(M2,2,nz,NumG)
  B3_14 = zeros(M2,2,nz,NumG)
  B1p_12 = zeros(2,nz,NumG)
  C23_2 = zeros(2,M2,nz,NumG)
  C14_3 = zeros(2,M2,nz,NumG)
  luSA = Array{LU}(undef,nz,NumG)
  SchurBand = Array{BandedMatrix}(undef,NumG)
  for iG = 1 : NumG
    SchurBand[iG] = BandedMatrix(FillArrays.Zeros(nz*4,nz*4), (3,3))  
  end  

  return Interior(
    M,
    nz,
    fac,
    FacGrav,
    A13,
    A23,
    A32,
    B1m_34,
    B1_1,
    B1_23,
    B1_4,
    B2_23,
    B3_14,
    B1p_12,
    C23_2,
    C14_3,
    luSA,
    SchurBand,
  )
end  

function FillInterior!(AI,U,DG,dz,fac,Phys,Param)
  M2 = AI.M - 2
  M1 = M - 1
  M23= M2 * 2 + M
  nz = AI.nz
  AI.fac = fac
  AI.fac = fac
  AI.FacGrav = 0.5 * Phys.Grav
  AI.FacGrav = 0.5 * Phys.Grav
  SA = zeros(M2,M2)
  SAD = zeros(M2,M2)
  DW = DG.DWZ
  wZ = DG.wZ
  wZ = DG.wZ
  ThPos = 5
  RhoPos = 1
  Th = zeros(M)
  dpdRhoTh = zeros(M)
  cS = Param.cS
  DoF  = size(U,3)
  S = zeros(2,2)

  for ID = 1 : DoF
    sh = 1
    @. AI.SchurBand[1].data = 0.0
    @. AI.SchurBand[1].data[4,:] = fac
    for iz = 1 : nz
      Th = U[:,iz,ID,ThPos]./U[:,iz,ID,RhoPos]
      dpdRhoTh = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[:,iz,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
      @views @. AI.A13[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[:,2:M1]
      @views AI.A23[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[2:M1,2:M1] * diagm(Th[2:M-1])
      @views AI.A32[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[2:M1,2:M1] * diagm(dpdRhoTh[2:M-1])
      @views SAD = fac * I - (0.5 * Phys.Grav / fac) * AI.A13[2:M-1,:,iz,ID] - 
        (1.0 / fac) *AI.A32[:,:,iz,ID] * AI.A23[:,:,iz,ID]
      AI.luSA[iz,ID] = LinearAlgebra.lu(SAD)

      if iz > 1
        AI.B1m_34[1,iz,ID] = -1.0 / (wZ[1] * dz[iz-1,ID])

        dpdRhoThM = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[M,iz-1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
        AI.B1m_34[2,iz,ID] = dpdRhoThM / (cS * wZ[1] * dz[iz-1,ID])
      end  
      @views @. AI.B1_23[:,1,iz,ID] = 2.0 * DW[:,1] / dz[iz,ID]
      @views @. AI.B1_23[:,2,iz,ID] = 2.0 *  DW[:,M] / dz[iz,ID]
      if iz == 1
        AI.B1_1[iz,ID] = 0
      else
        AI.B1_23[1,1,iz,ID] = 0.0
        AI.B1_1[iz,ID] = -dpdRhoTh[1] / (cS * wZ[1] * dz[iz,ID])
      end  
      if iz == nz
        AI.B1_4[iz,ID] = 0
      else
        AI.B1_23[M,2,iz,ID] = 0.0
        AI.B1_4[iz,ID] = -dpdRhoTh[M] / (cS * wZ[1] * dz[iz,ID])
      end  
      if iz < nz
        dpdRhoThP = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[1,iz+1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
        AI.B1p_12[1,iz,ID] = dpdRhoThP / (cS * wZ[1] * dz[iz+1])
        AI.B1p_12[2,iz,ID] = 1.0 / (wZ[1] * dz[iz+1])
      end  

      @views @. AI.B2_23[:,1,iz,ID] = 2.0 * DW[2:M1,1] * Th[1] / dz[iz,ID]
      @views @. AI.B2_23[:,2,iz,ID] = 2.0 * DW[2:M1,M] * Th[M] / dz[iz,ID]
      @views @. AI.B3_14[:,1,iz,ID] = 2.0 * DW[2:M1,1] * dpdRhoTh[1] / dz[iz,ID]
      @views @. AI.B3_14[:,2,iz,ID] = 2.0 * DW[2:M1,M] * dpdRhoTh[M] / dz[iz,ID]

      @views @. AI.C23_2[1,:,iz,ID] = 2.0 * DW[1,2:M1] * dpdRhoTh[2:M1] / dz[iz,ID]
      @views @. AI.C23_2[2,:,iz,ID] = 2.0 * DW[M,2:M1] * dpdRhoTh[2:M1] / dz[iz,ID]
      @views @. AI.C14_3[1,:,iz,ID] = 2.0 * DW[1,2:M1] * Th[2:M1] / dz[iz,ID]
      @views @. AI.C14_3[2,:,iz,ID] = 2.0 * DW[M,2:M1] * Th[2:M1] / dz[iz,ID]

      if iz > 1 
        ThM = U[M,iz-1,1,ThPos] / U[M,iz-1,1,RhoPos]  
        Th1 = U[1,iz-1,1,ThPos] / U[1,iz-1,1,RhoPos]  
        dpdRhoThM = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[M,iz-1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
        S[1,1] = ThM * dpdRhoThM / dz[iz-1,ID] / cS / wZ[1]
        S[2,1] = -Th[1] * dpdRhoThM / dz[iz-1,ID] / cS / wZ[1]
        S[1,2] = -ThM * dpdRhoTh[1] / dz[iz,ID] / cS / wZ[1]
        S[2,2] = Th[1] * dpdRhoTh[1] / dz[iz,ID] / cS / wZ[1]
        AI.SchurBand[ID][sh-1:sh,sh-1:sh] .+= S
        S[1,1] = cS / dz[iz-1,ID] / wZ[1]
        S[2,1] = -cS / dz[iz-1,ID] / wZ[1]
        S[1,2] = -cS / dz[iz,ID] / wZ[1]
        S[2,2] = cS / dz[iz,ID] / wZ[1]
        AI.SchurBand[ID][sh-2,sh-2] += S[1,1]
        AI.SchurBand[ID][sh-2,sh+1] += S[1,2]
        AI.SchurBand[ID][sh+1,sh-2] += S[2,1]
        AI.SchurBand[ID][sh+1,sh+1] += S[2,2]
        AI.SchurBand[ID][sh,sh-2] += -ThM / wZ[1] / dz[iz,ID]
        AI.SchurBand[ID][sh+1,sh-1] += -dpdRhoThM / wZ[1] / dz[iz,ID]
        AI.SchurBand[ID][sh-1,sh+1] += Th[1] / wZ[1] / dz[iz-1,ID]
        AI.SchurBand[ID][sh-2,sh] += dpdRhoTh[1] / wZ[1] / dz[iz-1,ID]
      end  
      if iz == 1
        AI.SchurBand[ID][sh,sh+1] += 2.0 * DW[1,1] * Th[1] / dz[iz,ID]
        AI.SchurBand[ID][sh+1,sh] += -2.0 * DW[1,1] * dpdRhoTh[1] / dz[iz,ID]
        AI.SchurBand[ID][sh+1,sh+1] += 2.0 * cS / dz[iz,ID] / wZ[1]  
      end  
      AI.SchurBand[ID][sh,sh+2] += 2.0 * DW[1,M] * Th[M] / dz[iz,ID]
      AI.SchurBand[ID][sh+1,sh+3] += 2.0 * DW[1,M] * dpdRhoTh[M] / dz[iz,ID]
      AI.SchurBand[ID][sh+3,sh+1] += 2.0 * DW[M,1] * Th[1] / dz[iz,ID]
      AI.SchurBand[ID][sh+2,sh] += 2.0 * DW[M,1] * dpdRhoTh[1] / dz[iz,ID]
      if iz == nz
        AI.SchurBand[ID][sh+3,sh+2] += 2.0 * DW[M,M] * Th[M] / dz[iz,ID]
        AI.SchurBand[ID][sh+2,sh+3] += -2.0 * DW[M,M] * dpdRhoTh[M] / dz[iz,ID]
        AI.SchurBand[ID][sh+2,sh+2] += 2.0 * cS / dz[iz,ID] / wZ[1]  
      end  
      sh += 4
    end
  end
end

function SchurBoundary!(AI)
  M2 = AI.M - 2
  invfac = 1.0 / AI.fac
  FacGrav = AI.FacGrav
  r1 = zeros(M,2)
  r2 = zeros(M2,2)
  r3 = zeros(M2,2)
  s = zeros(4,2)
  r11 = zeros(1,2)
  r1M = zeros(1,2)
  A22B = AI.SchurBand[1]
  DoF = size(AI.A13,4)
  for ID = 1 : DoF
    for iz = 1 : nz
      sh = (iz - 1) * 4  
#     Column 1 and 4
      @. r3 = AI.B3_14[:,:,iz,ID]   

      ldiv!(AI.luSA[iz,ID],r3)

      r11 = -invfac * (AI.A13[1:1,:,iz,ID] * r3[:,:])
      r1M = -invfac * (AI.A13[M:M,:,iz,ID] * r3[:,:])
      r11[1,1] += invfac * AI.B1_1[iz,ID]
      r1M[1,2] += invfac * AI.B1_4[iz,ID]
      @views r2 = -invfac * (AI.A23[:,:,iz,ID] * r3)

      @views s[[1,4],:] = -AI.C14_3[:,:,iz,ID] * r3
      @views s[[2,3],:] = -AI.C23_2[:,:,iz,ID] * r2
      s[2,1] = s[2,1] - FacGrav * r11[1,1]
      s[3,1] = s[3,1] - FacGrav * r1M[1,1]
      s[2,2] = s[2,2] - FacGrav * r11[1,2]
      s[3,2] = s[3,2] - FacGrav * r1M[1,2]
      A22B[sh + 1:sh + 4,[sh + 1,sh + 4]] .+= s

#     Column 2 and 3 
      @. r1 = AI.B1_23[:,:,iz,ID]
      @. r2 = AI.B2_23[:,:,iz,ID]
      @views r3 = -invfac * (AI.A32[:,:,iz,ID] * r2 + FacGrav * r1[2:M-1,:]) 

      ldiv!(AI.luSA[iz,ID],r3)

      r11 = invfac * (r1[1:1,:] - AI.A13[1:1,:,iz,ID] * r3[:,:])
      r1M = invfac * (r1[M:M,:] - AI.A13[M:M,:,iz,ID] * r3[:,:])
      @views r2 = invfac * (r2 - AI.A23[:,:,iz,ID] * r3)

      @views s[[1,4],:] = -AI.C14_3[:,:,iz,ID] * r3
      @views s[[2,3],:] = -AI.C23_2[:,:,iz,ID] * r2
      s[2,1] = s[2,1] - FacGrav * r11[1,1]
      s[3,1] = s[3,1] - FacGrav * r1M[1,1]
      s[2,2] = s[2,2] - FacGrav * r11[1,2]
      s[3,2] = s[3,2] - FacGrav * r1M[1,2]
      A22B[sh + 1:sh + 4,[sh + 2,sh + 3]] .+= s
      if iz > 1
        #Column -2   
        A22B[sh+2,sh-1] -= FacGrav * invfac * AI.B1m_34[1,iz,ID]
        # Column -1
        A22B[sh+2,sh] -= FacGrav * invfac * AI.B1m_34[2,iz,ID]
      end
      if iz < nz 
#       Column +1  
        A22B[sh+3,sh+5] -= FacGrav * invfac * AI.B1p_12[1,iz,ID]
#       Column +2  
        A22B[sh+3,sh+6] -= FacGrav * invfac * AI.B1p_12[2,iz,ID]
      end
    end    
    luBand!(A22B)
  end
end

function ldivBlockAF(invfac,FacGrav,A13,A23,A32,luSA,r1,r2,r3,r11,r1M)
   
  @views r3 .+= -invfac * (A32 * r2 + FacGrav * r1[2:end-1,:])

  ldiv!(luSA,r3)

  @views r11 .= invfac * (r1[1:1] - A13[1:1,:] * r3)
  @views r1M .= invfac * (r1[end:end] - A13[end:end,:] * r3)
  r2 .= invfac * (r2 - A23 * r3)
end    

function ldivBlockAB(invfac,FacGrav,A13,A23,A32,luSA,r1,r2,r3)

  @views r3 .+= -invfac * (A32 * r2 + FacGrav * r1[2:end-1,:])

  ldiv!(luSA,r3)

  r1 .= invfac * (r1 - A13 * r3)
  r2 .= invfac * (r2 - A23 * r3)
end

function ldivVertical!(AI,b)
  M2 = AI.M - 2
  M1 = AI.M - 1
  invfac = 1.0 / AI.fac
  FacGrav = AI.FacGrav
  r1 = zeros(M)
  r2 = zeros(M2)
  r3 = zeros(M2)
  rs = zeros(4*nz)
  s = zeros(4)
  r11 = zeros(1)
  r1M = zeros(1)
  A22B = AI.SchurBand[1]
  RhoPos = 1
  ThPos = 5
  wPos = 4
  DoF = size(b,3)
  for ID = 1 : DoF
#   Forward substitution  
    for iz = 1 : nz
      sh = (iz - 1) * 4  
      rs[sh + 1] = b[1,iz,ID,ThPos]
      rs[sh + 2] = b[1,iz,ID,wPos]
      rs[sh + 3] = b[M,iz,ID,wPos]
      rs[sh + 4] = b[M,iz,ID,ThPos]
      @views @. r1 = b[:,iz,ID,RhoPos]
      @views @. r2 = b[2:M1,iz,ID,ThPos]
      @views @. r3 = b[2:M1,iz,ID,wPos]

      @views ldivBlockAF(invfac,FacGrav,AI.A13[:,:,iz,ID],AI.A23[:,:,iz,ID],AI.A32[:,:,iz,ID],
        AI.luSA[iz,ID],r1,r2,r3,r11,r1M)

      @views s[[1,4]] = -AI.C14_3[:,:,iz,ID] * r3
      @views s[[2,3]] = -AI.C23_2[:,:,iz,ID] * r2
      s[2] = s[2] - FacGrav * r11[1]
      s[3] = s[3] - FacGrav * r1M[1]
      @views rs[sh + 1:sh + 4] .+= s
    end    
    ldivBand!(A22B,rs)
#   Back substitution  
    for iz = 1 : nz
      sh = (iz - 1) * 4
      b[1,iz,ID,ThPos] = rs[sh + 1]
      b[1,iz,ID,wPos] = rs[sh + 2]
      b[M,iz,ID,wPos] = rs[sh + 3]
      b[M,iz,ID,ThPos] = rs[sh + 4]
      @views @. r1 = b[:,iz,ID,RhoPos]
      @views @. r2 = b[2:M1,iz,ID,ThPos]
      @views @. r3 = b[2:M1,iz,ID,wPos]
      @views rsC = rs[sh + 1 : sh + 4] 
      r1[1] -= AI.B1_1[iz,ID] * rsC[1]
      r1[M] -= AI.B1_4[iz,ID] * rsC[4]
      @views r1 .-= AI.B1_23[:,:,iz,ID] * rsC[2:3]
      @views r2 .-= AI.B2_23[:,:,iz,ID] * rsC[2:3]
      @views r3 .-= AI.B3_14[:,:,iz,ID] * rsC[[1,4]]
      if iz > 1
        @views rsM = rs[sh - 1 : sh] 
        r1[1] -= AI.B1m_34[1,iz,ID] * rsM[1] + AI.B1m_34[2,iz,ID] * rsM[2]
      end
      if iz < nz
        @views rsP = rs[sh + 5 : sh + 6] 
        r1[M] -= AI.B1p_12[1,iz,ID] * rsP[1] + AI.B1p_12[2,iz,ID] * rsP[2]
      end
      @views ldivBlockAB(invfac,FacGrav,AI.A13[:,:,iz,ID],AI.A23[:,:,iz,ID],AI.A32[:,:,iz,ID],
        AI.luSA[iz,ID],r1,r2,r3)
      @views @. b[:,iz,ID,RhoPos] = r1
      @views @. b[2:M1,iz,ID,ThPos] = r2
      @views @. b[2:M1,iz,ID,wPos] = r3
    end  
  end
end

function Permutation(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  for iz = 1 : nz
    for iv = 1 : 1
      for k = 1 : M 
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
    for iv = 2 : 3
      for k = 2 : M - 1
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
  end
  ivw = 3
  ivTh = 2
  for iz = 1 : nz
      ii += 1
      p[ii] = 1 + (iz-1) * M  + (ivTh - 1) * N
      ii += 1
      p[ii] = 1 + (iz-1) * M  + (ivw - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivw - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
  end
  return p
end  
function PermutationAndres(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  ivRho = 1
  ivw = 3
  ivTh = 2
  for iz = 1 : nz
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivRho - 1) * N
    end
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivTh - 1) * N
    end
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivw - 1) * N
    end
  end
  for iz = 1 : nz
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivRho - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivw - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
  end
  return p
end

function PermutationAndres1(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  ivRho = 1
  ivw = 3
  ivTh = 2
  for iz = 1 : nz
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivRho - 1) * N
    end
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivTh - 1) * N
    end
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivw - 1) * N
    end
  end
  for iz = 1 : nz
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivRho - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivw - 1) * N
  end
  return p
end

function PermutationAndres2(M,nz)
#Permutation
  N = M * nz
  p = zeros(Int,3*N)
  ii = 0
  ivRho = 1
  ivw = 3
  ivTh = 2
  for iz = 1 : nz
    for k = 1 : M 
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivRho - 1) * N
    end
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivTh - 1) * N
    end
    for k = 1 : M - 1
      ii += 1
      p[ii] = k + (iz - 1) * M + (ivw - 1) * N
    end
  end
  for iz = 1 : nz
#     ii += 1
#     p[ii] = M + (iz - 1) * M + (ivRho - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
      ii += 1
      p[ii] = M + (iz - 1) * M + (ivw - 1) * N
  end
  return p
end


function Recover(JacS,M,nz,Physi,fac)
  M2 = M - 2
  m11 = M + 2 * (M - 2)
  A31 = [spzeros(M2,1) Phys.Grav*sparse(I,M2,M2) spzeros(M2,1) ]
  JacSP11 = spzeros(0,0)
  JacSP12 = spzeros(0,0)
  JacSP21 = spzeros(0,nz*m11)
  for iz = 1 : nz
    A13 = JacS.A13[:,:,iz,1]
    A23 = JacS.A23[:,:,iz,1]
    A32 = JacS.A32[:,:,iz,1]
    AS11=[fac*sparse(I,M,M)  spzeros(M,M2) A13
       spzeros(M2,M) fac*sparse(I,M2,M2) A23
       A31 A32 fac*sparse(I,M2,M2)]
    A31F = collect(A31)   
    SA = I(M2)*fac - (1.0 / fac) * A31F * A13 - (1.0 / fac) *A32 * A23
    DGSEMNeu.LUFull!(SA)
    @show sum(abs.(SA - JacS.SA[:,:,iz,1]))
    JacSP11 = [JacSP11 spzeros(m11*(iz-1),m11)
            spzeros(m11,m11*(iz-1)) AS11]
    B1_1 = JacS.B1_1[iz,1]
    B1_23 = JacS.B1_23[:,:,iz,1]
    B1_4 = JacS.B1_4[iz,1]
    B1m_34 = JacS.B1m_34[:,iz,1]'
    B1p_12 = JacS.B1p_12[:,iz,1]'
    B2_23 = JacS.B2_23[:,:,iz,1]
    B3_14  = JacS.B3_14[:,:,iz,1]
    BS12_l = [spzeros(m11,2) [B1m_34
                            spzeros(m11-1,2)]]
    BS12_r = [[spzeros(M-1,2)
               B1p_12
               spzeros(2*M2,2)] spzeros(m11,2)]
    BS12_c1 = [B1_1
               spzeros(M-1+M2,1)
               B3_14[:,1]] 
    BS12_c23 = [B1_23
                B2_23
                spzeros(M2,2)]
    BS12_c4 = [spzeros(M-1,1)
               B1_4
               spzeros(M2,1)
               B3_14[:,2]]
    BS12_c = [BS12_c1 BS12_c23 BS12_c4]           
    if iz == 1
      JacSP12 = [BS12_c BS12_r spzeros(m11,(nz-2)*4)]
    elseif iz == nz
      JacSP12 = [JacSP12
                 spzeros(m11,(nz-2)*4) BS12_l BS12_c]
    else             
      JacSP12 = [JacSP12
                 spzeros(m11,(iz-2)*4) BS12_l BS12_c BS12_r spzeros(m11,(nz-1-iz)*4)]
    end             
    C14_3 = JacS.C14_3[:,:,iz,1]
    C23_2 = JacS.C23_2[:,:,iz,1]
    CS11 = [spzeros(1,M)
            [Phys.Grav spzeros(1,M-1)]
            [spzeros(1,M-1) Phys.Grav ]
            spzeros(1,M)]
    CS12 = [spzeros(1,M2)
            C23_2
            spzeros(1,M2)]
    CS13 = [C14_3[1,:]'
            spzeros(2,M2)
            C14_3[2,:]']
    JacSP21 = [JacSP21
               [spzeros(4,(iz-1)*m11) CS11 CS12 CS13 spzeros(4,(nz-iz)*m11)]]

  end
  return JacSP11,JacSP12,JacSP21
end


backend = CPU()
FTB = Float64
Parallel = true

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

Problem = "WarmBubble2DXCart"
Param = Examples.Parameters(FTB,Problem)
# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

# ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Grid
Boundary = Grids.Boundary()
Boundary.WE = "Period"
Boundary.SN = "Period"
Boundary.BT = "FreeSlip"
nx = 2
ny = 2
Lx = 20000.0
Ly = 2000.0
H = 10000.0
x0 = 0.0
y0 = 0.0
nz = 10
OrdPoly = 4
OrdPolyZ = 6
M = OrdPolyZ + 1
N = M * nz
Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)
Grid.AdaptGrid = Grids.AdaptGrid(FTB,"Sleve",FTB(H))

OrdPrint = 4
OrdPrintZ = 4
DGMethod = ""
TopoS = ""
Topography=(TopoS=TopoS,
              )
TopoProfile = Examples.Flat()()
(DG, Metric, Exchange, Global) = DyCore.InitCartDG(backend,FTB,OrdPoly,OrdPolyZ,DGMethod,
    OrdPrint,OrdPrintZ,H,Topography,Model,
    Phys,TopoProfile,CellToProc,Grid,ParallelCom)
#DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)


fac = 0.5
U = ones(M,nz,DG.NumG,5)
Aux = rand(M,nz,DG.NumG,4)
@. Aux += 50.0
H = 10000.0
dzLoc = H / nz
dz = ones(nz,DG.NumG) * dzLoc
Metric.dz = dz

VelForm = Examples.VelocityC()
Examples.InitialProfile!(backend,FTB,Model,Problem,Param,Phys,VelForm)
Profile = Model.InitialProfile
time = 0.0
for iz = 1 : nz
  z0 = (iz - 1) * dzLoc  
  z1 = iz * dzLoc  
  for k = 1 : M
    zLoc = 0.5 * ((1.0 - DG.xwZ[k]) * z0 + (1.0 + DG.xwZ[k]) * z1)  
    xS = SVector{3}(0.0,0.0,zLoc)
    RhoP,_,_,_,ThP= Profile(xS,time)
    @views @. U[k,iz,:,1] = RhoP
    @views @. U[k,iz,:,5] = RhoP * ThP
  end
end

Model.GPAuxPos = 2
Model.GeoPotential = Sources.GeoPotentialDeep()(Model.GPAuxPos,Grid.Form)
DGSEMNeu.GeoPot(Aux,DG,Metric,Exchange,Global)
Geo = Aux[:,:,:,Model.GPAuxPos]

# Version Marco
dt = 2.0
JCacheA = DGSEMNeu.JacobianCacheMarcoSplitNS(backend,FTB,M,nz,DG)
DGSEMNeu.precompute_gravity!(Geo,Metric.dz,JCacheA)
DGSEMNeu.FillJacDGVertKernel!(JCacheA, U, Metric.dz, fac)
DGSEMNeu.SchurBoundaryKernel!(JCacheA, Metric.dz, Phys)
DGSEMNeu.luBandkernel!(JCacheA)
b = rand(M,nz,DG.NumG,5)
c = zeros(M,nz,DG.NumG,5)
d = zeros(M,nz,DG.NumG,5)
@. b = b + 5.0
@. c = b
@. d = b

DGSEMNeu.solve_jacobian!(b, JCacheA, Metric)

# Version Oswald
dt = 2.0
JCacheO = DGSEMNeu.JacDGVert(backend,FTB,M,nz,DG)
DGSEMNeu.FillJacDGVert!(JCacheO, U, DG, Metric.dz, fac, Phys)
DGSEMNeu.SchurBoundary!(JCacheO)
DGSEMNeu.Solve!(JCacheO,c)


Gblk, dGeo = DGSEMNeu.DerivativeGeo(Geo,DG,dz)
dSdS,dSdM,dMdS,dMdM = DGSEMNeu.InitJacDG(DG,nz,Param)

fac = 0.5
JacLU, Jac = DGSEMNeu.JacDGTNeu(U,Aux,DG,fac,dSdS,dSdM,dMdS,dMdM,dGeo,dz,Phys)


dF = zeros(M*nz*3)
idF = 0
@inbounds for iv in [1,5,4]
  @inbounds for iz = 1 : nz
    @inbounds for i = 1 : M
      global idF += 1
      dF[idF] = d[i,iz,1,iv]
    end
  end
end
ldiv!(JacLU[1],dF)
idF = 0
@inbounds for iv in [1,5,4]
  @inbounds for iz = 1 : nz
    @inbounds for i = 1 : M
      global idF += 1
      d[i,iz,1,iv] = dF[idF] 
    end
  end
end
stop

DGSEMNeu.FillJacDGVert!(JacS,U,DG,dz,fac,Phys)

#DGSEMNeu.SchurBoundary!(JacS)

p = PermutationAndres1(M,nz)
JacP = Jac[p,p]
n22 = 3*nz
n11 = 3*M*nz - n22
JacP11 = JacP[1:n11,1:n11]
JacP12 = JacP[1:n11,n11+1:end]
JacP21 = JacP[n11+1:end,1:n11]
JacP22 = JacP[n11+1:end,n11+1:end]
JacP22B = BandedMatrix(JacP22)

JacP11F = collect(JacP11)
JacP12F = collect(JacP12)
JacP21F = collect(JacP21)
JacP22F = collect(JacP22)

SchurJacPF = JacP22F - JacP21F * (JacP11F \ JacP12F)
SchurJacP = sparse(SchurJacPF)
SchurJacPB = BandedMatrix(SchurJacP)
stop



#BP11,BP12,BP21,BP22= Recover(JacS,M,nz,Phys,fac)

b = ones(size(U))
@. b[:,:,:,5] *= 2
@. b[:,:,:,4] *= 3

bF = zeros(M*nz*3)
ibF = 0
@inbounds for iv in [1,4,5]
  @inbounds for iz = 1 : nz
    @inbounds for i = 1 : M
      global ibF += 1
      bF[ibF] = b[i,iz,1,iv]
    end
  end
end
JacVLU, A = DGSEMNeu.JacDG(U,Aux,DG,fac,dSdS,dSdM,dMdS,dMdM,dz,Phys)
#JacVLU = DGSEMNeu.JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,dz,Phys)
ldiv!(JacVLU[1],bF)
DGSEMNeu.ldivVertical!(JacS,b)

stop
bR = zeros(M*nz*3)
bRP = zeros(M*nz*3)
b1R = reshape(b1,M*nz*3)
@. bR = b1R
bRP .= bR[p]
ldiv!(JacLU[1],bR)
DGSEMNeu.ldivVertical!(JacS,b)
@views @. b1[:,:,1] = b[:,:,1,1]
@views @. b1[:,:,2] = b[:,:,1,5]
@views @. b1[:,:,3] = b[:,:,1,4]
xP = JacP \ bRP

xP1 = JacP11 \ bRP[1:n11]
rP2 = bRP[n11+1:end] - JacP21 * xP1


aaa=3
