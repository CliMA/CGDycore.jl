import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGVertical, GPU, DyCore, DGSEM
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

Problem = "HillAgnesiXCart"
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
Lx = 1000.0
Ly = 1000.0
x0 = 0.0
y0 = 0.0
nz = 5
OrdPoly = 4
OrdPolyZ = 7
M = OrdPolyZ + 1
N = M * nz
Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)

OrdPrint = 4
OrdPrintZ = 4
DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)


fac = 0.5
U = ones(M,nz,DG.NumG,5)
H = 10000.0
dzLoc = H / nz
dz = ones(nz,DG.NumG) * dzLoc

Profile = Examples.StratifiedExample()(Param,Phys)
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

AI = Interior(M,nz,size(U,3))
FillInterior!(AI,U,DG,dz,fac,Phys,Param)
SchurBoundary!(AI)

b = ones(size(U))
@. b[:,:,:,5] *= 2
@. b[:,:,:,4] *= 3
ldivVertical!(AI,b)

