mutable struct JacDGVert 
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
  SA::Array{Float64, 4}
  SchurBand::Array{Float64, 3}
end  

function LUFull!(A)

  n = size(A,1)
  @inbounds for j = 1 : n - 1
    @inbounds for i = j + 1 : n 
      A[i,j] /= A[j,j]  
      @inbounds for k = j + 1 : n
        A[i,k] -= A[i,j] * A[j,k]
      end  
    end  
  end  
end

function ldivFull2!(A,b)

# Forward loop
  n = size(A,1)
  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : n
      @views @.  b[i,:] -= A[i,k] * b[k,:]
    end
  end
#  Backward loop
  @inbounds for k = n : -1 : 1
    @views @.  b[k,:] /= A[k,k]
    @inbounds for i = 1 : k - 1
      @views @.  b[i,:] -= A[i,k] * b[k,:]
    end
  end
end

function ldivFull!(A,b)

# Forward loop
  n = size(A,1)
  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : n
      @views b[i] -= A[i,k] * b[k]
    end
  end
#  Backward loop
  @inbounds for k = n : -1 : 1
    @views b[k] /= A[k,k]
    @inbounds for i = 1 : k - 1
      @views b[i] -= A[i,k] * b[k]
    end
  end
end

function ldivBand!(A,b,l,u)

# Ly = b
  A = A
  n = size(A,2)
  up1 = u + 1 
  uplp1 = up1 + l

  @inbounds for k = 1 : n - 1
    @inbounds for i = k + 1 : min(k+l,n)
      b[i] -= A[i - k + up1,k] * b[k]  
    end    
  end
# Ux = y
  @inbounds for k = n : -1 : 1
    b[k] /= A[up1,k]
    @inbounds for i = max(k - u, 1) : k - 1
      b[i] -= A[up1 + i - k,k] * b[k]  
    end    
  end
end

function luBand!(A,l,u)
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
  n = size(A,2)
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

function JacDGVert(M,nz,NumG)
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
  SA = zeros(M2,M2,nz,NumG)
  SchurBand = zeros(7,4*nz,NumG)

  return JacDGVert(
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
    SA,
    SchurBand,
  )
end  

function FillJacDGVert!(JacVert,U,DG,dz,fac,Phys,Param)
  backend = get_backend(U)
  FTB = eltype(U)
  M = JacVert.M
  M2 = JacVert.M - 2
  M1 = M - 1
  M23= M2 * 2 + M
  nz = JacVert.nz
  JacVert.fac = fac
  JacVert.fac = fac
  JacVert.FacGrav = 0.5 * Phys.Grav
  JacVert.FacGrav = 0.5 * Phys.Grav
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

  @inbounds for ID = 1 : DoF
    sh = 1
    @views @. JacVert.SchurBand[:,:,ID] = 0.0
    @views @. JacVert.SchurBand[4,:,ID] = fac
    @time @inbounds for iz = 1 : nz
      @views @. Th = U[:,iz,ID,ThPos]/U[:,iz,ID,RhoPos]
      @views @. dpdRhoTh = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
        (Phys.Rd * U[:,iz,ID,ThPos] / Phys.p0)^(Phys.kappa / (1.0 - Phys.kappa))  
      @views @. JacVert.A13[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[:,2:M1]
      #@views JacVert.A23[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[2:M1,2:M1] * diagm(Th[2:M-1])
      #@views JacVert.A32[:,:,iz,ID] = (2 / dz[iz,ID]) * DW[2:M1,2:M1] * diagm(dpdRhoTh[2:M-1])
      @inbounds for j = 2 : M1
        @views @. JacVert.A23[:,j-1,iz,ID] = 2 / dz[iz,ID] * DW[j,2:M1] * Th[j]  
        @views @. JacVert.A32[:,j-1,iz,ID] = 2 / dz[iz,ID] * DW[j,2:M1] * dpdRhoTh[j]  
      end  
      #@views JacVert.SA[:,:,iz,ID] .= fac * I - (0.5 * Phys.Grav / fac) * JacVert.A13[2:M-1,:,iz,ID] - 
      #  (1.0 / fac) *JacVert.A32[:,:,iz,ID] * JacVert.A23[:,:,iz,ID]
      @inbounds for i = 1 : M2
        @inbounds for j = 1 : M2
          JacVert.SA[i,j,iz,ID] = 0.0
          @inbounds for k = 1 : M2
            JacVert.SA[i,j,iz,ID] += JacVert.A32[i,k,iz,ID] * JacVert.A23[k,j,iz,ID]
          end
          if i == j
            JacVert.SA[i,j,iz,ID] = fac - (0.5 * Phys.Grav / fac) * JacVert.A13[i+1,j,iz,ID] -
              (1.0 / fac) * JacVert.SA[i,j,iz,ID]
          else
            JacVert.SA[i,j,iz,ID] = -(0.5 * Phys.Grav / fac) * JacVert.A13[i+1,j,iz,ID] -
              (1.0 / fac) * JacVert.SA[i,j,iz,ID]
          end    
        end
      end  
      @views LUFull!(JacVert.SA[:,:,iz,ID])  

      if iz > 1
        JacVert.B1m_34[1,iz,ID] = -1.0 / (wZ[1] * dz[iz-1,ID])

        dpdRhoThM = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[M,iz-1,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa))  
        JacVert.B1m_34[2,iz,ID] = dpdRhoThM / (cS * wZ[1] * dz[iz-1,ID])
      end  
      @views @. JacVert.B1_23[:,1,iz,ID] = 2.0 * DW[:,1] / dz[iz,ID]
      @views @. JacVert.B1_23[:,2,iz,ID] = 2.0 *  DW[:,M] / dz[iz,ID]
      if iz == 1
        JacVert.B1_1[iz,ID] = 0
      else
        JacVert.B1_23[1,1,iz,ID] = 0.0
        JacVert.B1_1[iz,ID] = -dpdRhoTh[1] / (cS * wZ[1] * dz[iz,ID])
      end  
      if iz == nz
        JacVert.B1_4[iz,ID] = 0
      else
        JacVert.B1_23[M,2,iz,ID] = 0.0
        JacVert.B1_4[iz,ID] = -dpdRhoTh[M] / (cS * wZ[1] * dz[iz,ID])
      end  
      if iz < nz
        dpdRhoThP = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[1,iz+1,ID,ThPos] / Phys.p0)^(Phys.kappa / (1.0 - Phys.kappa))  
        JacVert.B1p_12[1,iz,ID] = dpdRhoThP / (cS * wZ[1] * dz[iz+1,ID])
        JacVert.B1p_12[2,iz,ID] = 1.0 / (wZ[1] * dz[iz+1,ID])
      end  

      @views @. JacVert.B2_23[:,1,iz,ID] = 2.0 * DW[2:M1,1] * Th[1] / dz[iz,ID]
      @views @. JacVert.B2_23[:,2,iz,ID] = 2.0 * DW[2:M1,M] * Th[M] / dz[iz,ID]
      @views @. JacVert.B3_14[:,1,iz,ID] = 2.0 * DW[2:M1,1] * dpdRhoTh[1] / dz[iz,ID]
      @views @. JacVert.B3_14[:,2,iz,ID] = 2.0 * DW[2:M1,M] * dpdRhoTh[M] / dz[iz,ID]

      @views @. JacVert.C23_2[1,:,iz,ID] = 2.0 * DW[1,2:M1] * dpdRhoTh[2:M1] / dz[iz,ID]
      @views @. JacVert.C23_2[2,:,iz,ID] = 2.0 * DW[M,2:M1] * dpdRhoTh[2:M1] / dz[iz,ID]
      @views @. JacVert.C14_3[1,:,iz,ID] = 2.0 * DW[1,2:M1] * Th[2:M1] / dz[iz,ID]
      @views @. JacVert.C14_3[2,:,iz,ID] = 2.0 * DW[M,2:M1] * Th[2:M1] / dz[iz,ID]

      if iz > 1 
        ThM = U[M,iz-1,1,ThPos] / U[M,iz-1,1,RhoPos]  
        Th1 = U[1,iz-1,1,ThPos] / U[1,iz-1,1,RhoPos]  
        dpdRhoThM = FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
          (Phys.Rd * U[M,iz-1,ID,ThPos] / Phys.p0)^(Phys.kappa / (1.0 - Phys.kappa))  
        S[1,1] = ThM * dpdRhoThM / dz[iz-1,ID] / cS / wZ[1]
        S[2,1] = -Th[1] * dpdRhoThM / dz[iz-1,ID] / cS / wZ[1]
        S[1,2] = -ThM * dpdRhoTh[1] / dz[iz,ID] / cS / wZ[1]
        S[2,2] = Th[1] * dpdRhoTh[1] / dz[iz,ID] / cS / wZ[1]
        #JacVert.SchurBand[ID][sh-1:sh,sh-1:sh] .+= S
        #JacVert.SchurBand[ID][sh-1,sh-1] += S[1,1]
        JacVert.SchurBand[(sh-1) - (sh-1) + 4,sh-1,ID] += S[1,1]
        #JacVert.SchurBand[ID][sh  ,sh-1] += S[2,1]
        JacVert.SchurBand[sh - (sh-1) + 4,sh-1,ID] += S[2,1]
        #JacVert.SchurBand[ID][sh-1,sh  ] += S[1,2]
        JacVert.SchurBand[(sh-1) - sh + 4,sh,ID] += S[1,2]
        #JacVert.SchurBand[ID][sh  ,sh  ] += S[2,2]
        JacVert.SchurBand[sh - sh + 4,sh,ID] += S[2,2]
        S[1,1] = cS / dz[iz-1,ID] / wZ[1]
        S[2,1] = -cS / dz[iz-1,ID] / wZ[1]
        S[1,2] = -cS / dz[iz,ID] / wZ[1]
        S[2,2] = cS / dz[iz,ID] / wZ[1]
        #JacVert.SchurBand[ID][sh-2,sh-2] += S[1,1]
        JacVert.SchurBand[(sh-2) - (sh-2) + 4,sh-2,ID] += S[1,1]
        #JacVert.SchurBand[ID][sh-2,sh+1] += S[1,2]
        JacVert.SchurBand[(sh-2) - (sh+1) + 4,sh+1,ID] += S[1,2]
        #JacVert.SchurBand[ID][sh+1,sh-2] += S[2,1]
        JacVert.SchurBand[(sh+1) - (sh-2) + 4,sh-2,ID] += S[2,1]
        #JacVert.SchurBand[ID][sh+1,sh+1] += S[2,2]
        JacVert.SchurBand[(sh+1) - (sh+1) + 4,sh+1,ID] += S[2,2]
        #JacVert.SchurBand[ID][sh,sh-2] += -ThM / wZ[1] / dz[iz,ID]
        JacVert.SchurBand[sh - (sh-2) + 4,sh-2,ID] += -ThM / wZ[1] / dz[iz,ID]
        #JacVert.SchurBand[ID][sh+1,sh-1] += -dpdRhoThM / wZ[1] / dz[iz,ID]
        JacVert.SchurBand[(sh+1) - (sh-1) + 4,sh-1,ID] += -dpdRhoThM / wZ[1] / dz[iz,ID]
        #JacVert.SchurBand[ID][sh-1,sh+1] += Th[1] / wZ[1] / dz[iz-1,ID]
        JacVert.SchurBand[(sh-1) - (sh+1) + 4,sh+1,ID] += Th[1] / wZ[1] / dz[iz-1,ID]
        #JacVert.SchurBand[ID][sh-2,sh] += dpdRhoTh[1] / wZ[1] / dz[iz-1,ID]
        JacVert.SchurBand[(sh-2) - sh + 4,sh,ID] += dpdRhoTh[1] / wZ[1] / dz[iz-1,ID]
      end  
      if iz == 1
        #JacVert.SchurBand[ID][sh,sh+1] += 2.0 * DW[1,1] * Th[1] / dz[iz,ID]
        JacVert.SchurBand[sh - (sh+1) + 4,sh+1,ID] += 2.0 * DW[1,1] * Th[1] / dz[iz,ID]
        #JacVert.SchurBand[ID][sh+1,sh] += -2.0 * DW[1,1] * dpdRhoTh[1] / dz[iz,ID]
        JacVert.SchurBand[(sh+1) - sh + 4,sh,ID] += -2.0 * DW[1,1] * dpdRhoTh[1] / dz[iz,ID]
        #JacVert.SchurBand[ID][sh+1,sh+1] += 2.0 * cS / dz[iz,ID] / wZ[1]  
        JacVert.SchurBand[(sh+1) - (sh+1) + 4,sh+1,ID] += 2.0 * cS / dz[iz,ID] / wZ[1]
      end  
      #JacVert.SchurBand[ID][sh,sh+2] += 2.0 * DW[1,M] * Th[M] / dz[iz,ID]
      JacVert.SchurBand[sh - (sh+2) + 4,sh+2,ID] += 2.0 * DW[1,M] * Th[M] / dz[iz,ID]
      #JacVert.SchurBand[ID][sh+1,sh+3] += 2.0 * DW[1,M] * dpdRhoTh[M] / dz[iz,ID]
      JacVert.SchurBand[(sh+1) - (sh+3) + 4,sh+3,ID] += 2.0 * DW[1,M] * dpdRhoTh[M] / dz[iz,ID]
      #JacVert.SchurBand[ID][sh+3,sh+1] += 2.0 * DW[M,1] * Th[1] / dz[iz,ID]
      JacVert.SchurBand[(sh+3) - (sh+1) + 4,sh+1,ID] += 2.0 * DW[M,1] * Th[1] / dz[iz,ID]
      #JacVert.SchurBand[ID][sh+2,sh] += 2.0 * DW[M,1] * dpdRhoTh[1] / dz[iz,ID]
      JacVert.SchurBand[(sh+2) - sh + 4,sh,ID] += 2.0 * DW[M,1] * dpdRhoTh[1] / dz[iz,ID]
      if iz == nz
        #JacVert.SchurBand[ID][sh+3,sh+2] += 2.0 * DW[M,M] * Th[M] / dz[iz,ID]
        JacVert.SchurBand[(sh+3) - (sh+2) + 4,sh+2,ID] += 2.0 * DW[M,M] * Th[M] / dz[iz,ID]
        #JacVert.SchurBand[ID][sh+2,sh+3] += -2.0 * DW[M,M] * dpdRhoTh[M] / dz[iz,ID]
        JacVert.SchurBand[(sh+2) - (sh+3) + 4,sh+3,ID] += -2.0 * DW[M,M] * dpdRhoTh[M] / dz[iz,ID]
        #JacVert.SchurBand[ID][sh+2,sh+2] += 2.0 * cS / dz[iz,ID] / wZ[1]  
        JacVert.SchurBand[(sh+2) - (sh+2) + 4,sh+2,ID] += 2.0 * cS / dz[iz,ID] / wZ[1]
      end  
      sh += 4
    end
  end
  stop
end

function SchurBoundary!(JacVert)
  M = JacVert.M
  nz = JacVert.nz
  M2 = M - 2
  invfac = 1.0 / JacVert.fac
  FacGrav = JacVert.FacGrav
  r1 = zeros(M,2)
  r2 = zeros(M2,2)
  r3 = zeros(M2,2)
  s = zeros(4,2)
  r11 = zeros(2)
  r1M = zeros(2)
  DoF = size(JacVert.A13,4)
  @inbounds for ID = 1 : DoF
    @views A22B = JacVert.SchurBand[:,:,ID]
    @inbounds for iz = 1 : nz
      sh = (iz - 1) * 4  
#     Column 1 and 4
      @views @. r3 = JacVert.B3_14[:,:,iz,ID]   

      @views ldivFull2!(JacVert.SA[:,:,iz,ID],r3)

      #r11 = -invfac * (JacVert.A13[1:1,:,iz,ID] * r3[:,:])
      #r1M = -invfac * (JacVert.A13[M:M,:,iz,ID] * r3[:,:])
      #r11[1,1] += invfac * JacVert.B1_1[iz,ID]
      #r1M[1,2] += invfac * JacVert.B1_4[iz,ID]

      r11[1] = -JacVert.B1_1[iz,ID]
      r11[2] = 0.0
      r1M[1] = 0.0
      r1M[2] = -JacVert.B1_4[iz,ID]
      @inbounds for j = 1 : M2
        @views @. r11 += JacVert.A13[1,j,iz,ID] * r3[j,:] 
        @views @. r1M += JacVert.A13[M,j,iz,ID] * r3[j,:] 
      end    
      r11 .*= invfac
      r1M .*= invfac

      #@views r2 = -invfac * (JacVert.A23[:,:,iz,ID] * r3)
      @inbounds for i = 1 : M2
        @views @. r2[i,:] = 0.0 
        @inbounds for j = 1 : M2
          @views @. r2[i,:] += JacVert.A23[i,j,iz,ID] * r3[j,:]
        end  
        @views @. r2[i,:] *= -invfac
      end  

      #@views s[[1,4],:] = -JacVert.C14_3[:,:,iz,ID] * r3
      #@views s[[2,3],:] = -JacVert.C23_2[:,:,iz,ID] * r2
      @. s = 0.0
      @inbounds for j = 1 : M2
        @views @. s[1,:] += -JacVert.C14_3[1,j,iz,ID] * r3[j,:]
        @views @. s[4,:] += -JacVert.C14_3[2,j,iz,ID] * r3[j,:]
        @views @. s[2,:] += -JacVert.C23_2[1,j,iz,ID] * r2[j,:]
        @views @. s[3,:] += -JacVert.C23_2[2,j,iz,ID] * r2[j,:]
      end
      s[2,1] = s[2,1] - FacGrav * r11[1]
      s[3,1] = s[3,1] - FacGrav * r1M[1]
      s[2,2] = s[2,2] - FacGrav * r11[2]
      s[3,2] = s[3,2] - FacGrav * r1M[2]
      #@views A22B[sh + 1:sh + 4,[sh + 1,sh + 4]] .+= s
      @inbounds for i = 1 : 4
         A22B[i + sh - (sh +1) + 4,sh + 1] += s[i,1]
         A22B[i + sh - (sh +4) + 4,sh + 4] += s[i,2]
      end   

#     Column 2 and 3 
      @views @. r1 = JacVert.B1_23[:,:,iz,ID]
      @views @. r2 = JacVert.B2_23[:,:,iz,ID]
      #@views r3 = -invfac * (JacVert.A32[:,:,iz,ID] * r2 + FacGrav * r1[2:M-1,:]) 
      @inbounds for i = 1 : M2
        @views @. r3[i,:] = FacGrav * r1[i+1,:]
        @inbounds for j = 1 : M2
          @views @. r3[i,:] += JacVert.A32[i,j,iz,ID] * r2[j,:] 
        end  
        @views @. r3[i,:] *= -invfac
      end  

      @views ldivFull2!(JacVert.SA[:,:,iz,ID],r3)

      #r11 = invfac * (r1[1:1,:] - JacVert.A13[1:1,:,iz,ID] * r3[:,:])
      #r1M = invfac * (r1[M:M,:] - JacVert.A13[M:M,:,iz,ID] * r3[:,:])
      @views @. r11 = r1[1,:]
      @views @. r1M = r1[M,:]
      @inbounds for j = 1 : M2
        @views @. r11 += -JacVert.A13[1,j,iz,ID] * r3[j,:]
        @views @. r1M += -JacVert.A13[M,j,iz,ID] * r3[j,:]
      end  
      @views @. r11 *= invfac
      @views @. r1M *= invfac

      #@views r2 = invfac * (r2 - JacVert.A23[:,:,iz,ID] * r3)
      @inbounds for i = 1 : M2
        @inbounds for j = 1 : M2
          @views @. r2[i,:] += -JacVert.A23[i,j,iz,ID] * r3[j,:]
        end
        @views @. r2[i,:] *= invfac
      end  
        

      #@views s[[1,4],:] = -JacVert.C14_3[:,:,iz,ID] * r3
      #@views s[[2,3],:] = -JacVert.C23_2[:,:,iz,ID] * r2
      @. s = 0.0
      @inbounds for j = 1 : M2
        @views @. s[1,:] += - JacVert.C14_3[1,j,iz,ID] * r3[j,:]
        @views @. s[4,:] += - JacVert.C14_3[2,j,iz,ID] * r3[j,:]
        @views @. s[2,:] += - JacVert.C23_2[1,j,iz,ID] * r2[j,:]
        @views @. s[3,:] += - JacVert.C23_2[2,j,iz,ID] * r2[j,:]
      end  
        
      s[2,1] = s[2,1] - FacGrav * r11[1]
      s[3,1] = s[3,1] - FacGrav * r1M[1]
      s[2,2] = s[2,2] - FacGrav * r11[2]
      s[3,2] = s[3,2] - FacGrav * r1M[2]
      #@views A22B[sh + 1:sh + 4,[sh + 2,sh + 3]] .+= s
      @inbounds for i = 1 : 4
         A22B[i + sh - (sh +2) + 4,sh + 2] += s[i,1]
         A22B[i + sh - (sh +3) + 4,sh + 3] += s[i,2]
      end   
      if iz > 1
        #Column -2   
        #@time A22B[sh+2,sh-1] -= FacGrav * invfac * JacVert.B1m_34[1,iz,ID]
        A22B[sh+2-(sh-1)+4,sh-1] -= FacGrav * invfac * JacVert.B1m_34[1,iz,ID]
        # Column -1
        #A22B[sh+2,sh] -= FacGrav * invfac * JacVert.B1m_34[2,iz,ID]
        A22B[sh+2-sh+4,sh] -= FacGrav * invfac * JacVert.B1m_34[2,iz,ID]
      end
      if iz < nz 
#       Column +1  
        #A22B[sh+3,sh+5] -= FacGrav * invfac * JacVert.B1p_12[1,iz,ID]
        A22B[sh+3-(sh+5)+4,sh+5] -= FacGrav * invfac * JacVert.B1p_12[1,iz,ID]
#       Column +2  
        #A22B[sh+3,sh+6] -= FacGrav * invfac * JacVert.B1p_12[2,iz,ID]
        A22B[sh+3-(sh+6)+4,sh+6] -= FacGrav * invfac * JacVert.B1p_12[2,iz,ID]
      end
    end    
    luBand!(A22B,3,3)
  end
end

function ldivBlockAF(M,invfac,FacGrav,A13,A23,A32,SA,r1,r2,r3,r11,r1M)
   
  @inbounds for i = 1 : M - 2 
    r3[i] = FacGrav * r1[i+1]
    @inbounds for j = 1 : M - 2  
      r3[i] += A32[i,j] * r2[j]
    end
    r3[i] *= -invfac
  end  

  ldivFull!(SA,r3)

  r11[1,1] = r1[1]
  r1M[1,1] = r1[M]
  @inbounds for j = 1 : M - 2
    r11[1,1] += -A13[1,j] * r3[j]  
    r1M[1,1] += -A13[M,j] * r3[j]  
  end  
  r11[1,1] *= invfac
  r1M[1,1] *= invfac

  @inbounds for i = 1 : M - 2
    @inbounds for j = 1 : M - 2  
      r2[i] += -A23[i,j] * r3[j]  
    end
    r2[i] *= invfac
  end  
end    

function ldivBlockAB(M,invfac,FacGrav,A13,A23,A32,SA,r1,r2,r3)

  #@views r3 .+= -invfac * (A32 * r2 + FacGrav * r1[2:end-1])
  @inbounds for i = 1 : M - 2 
    r3[i] += -invfac * FacGrav * r1[i+1]
    @inbounds for j = 1 : M - 2  
      r3[i] += -invfac * A32[i,j] * r2[j]
    end
  end  

  ldivFull!(SA,r3)

  #r1 .= invfac * (r1 - A13 * r3)
  @inbounds for i = 1 : M
    @inbounds for j = 1 : M - 2  
      r1[i] += -A13[i,j] * r3[j]  
    end
    r1[i] *= invfac
  end  
  #r2 .= invfac * (r2 - A23 * r3)
  @inbounds for i = 1 : M - 2
    @inbounds for j = 1 : M - 2  
      r2[i] += -A23[i,j] * r3[j]  
    end
    r2[i] *= invfac
  end  
end

function ldivVertical!(JacVert,b)
  M = JacVert.M
  nz = JacVert.nz
  M2 = M - 2
  M1 = M - 1
  invfac = 1.0 / JacVert.fac
  FacGrav = JacVert.FacGrav
  r1 = zeros(M)
  r2 = zeros(M2)
  r3 = zeros(M2)
  rs = zeros(4*nz)
  s = zeros(4)
  r11 = zeros(1)
  r1M = zeros(1)
  RhoPos = 1
  ThPos = 5
  wPos = 4
  DoF = size(b,3)
  @inbounds for ID = 1 : DoF
    @views A22B = JacVert.SchurBand[:,:,ID]
#   Forward substitution  
    @inbounds for iz = 1 : nz
      sh = (iz - 1) * 4  
      rs[sh + 1] = b[1,iz,ID,ThPos]
      rs[sh + 2] = b[1,iz,ID,wPos]
      rs[sh + 3] = b[M,iz,ID,wPos]
      rs[sh + 4] = b[M,iz,ID,ThPos]
      @views @. r1 = b[:,iz,ID,RhoPos]
      @views @. r2 = b[2:M1,iz,ID,ThPos]
      @views @. r3 = b[2:M1,iz,ID,wPos]

      @views ldivBlockAF(M,invfac,FacGrav,JacVert.A13[:,:,iz,ID],JacVert.A23[:,:,iz,ID],JacVert.A32[:,:,iz,ID],
        JacVert.SA[:,:,iz,ID],r1,r2,r3,r11,r1M)

      @inbounds for j = 1 : M2
        rs[sh + 1] += -JacVert.C14_3[1,j,iz,ID] * r3[j]
        rs[sh + 4] += -JacVert.C14_3[2,j,iz,ID] * r3[j]
        rs[sh + 2] += -JacVert.C23_2[1,j,iz,ID] * r2[j]
        rs[sh + 3] += -JacVert.C23_2[2,j,iz,ID] * r2[j]
      end  
      rs[sh + 2] += -FacGrav * r11[1]
      rs[sh + 3] += -FacGrav * r1M[1]
    end    
    ldivBand!(A22B,rs,3,3)
#   Back substitution  
    @inbounds for iz = 1 : nz
      sh = (iz - 1) * 4
      b[1,iz,ID,ThPos] = rs[sh + 1]
      b[1,iz,ID,wPos] = rs[sh + 2]
      b[M,iz,ID,wPos] = rs[sh + 3]
      b[M,iz,ID,ThPos] = rs[sh + 4]
      @views @. r1 = b[:,iz,ID,RhoPos]
      @views @. r2 = b[2:M1,iz,ID,ThPos]
      @views @. r3 = b[2:M1,iz,ID,wPos]
      #@views rsC = rs[sh + 1 : sh + 4] 
      #r1[1] -= JacVert.B1_1[iz,ID] * rsC[1]
      #r1[M] -= JacVert.B1_4[iz,ID] * rsC[4]
      r1[1] -= JacVert.B1_1[iz,ID] * rs[sh + 1]
      r1[M] -= JacVert.B1_4[iz,ID] * rs[sh + 4]
      #@views r1 .-= JacVert.B1_23[:,:,iz,ID] * rsC[2:3]
      for j = 1 : M
        r1[j] -= JacVert.B1_23[j,1,iz,ID] * rs[sh + 2] + JacVert.B1_23[j,2,iz,ID] * rs[sh + 3] 
      end  
      #@views r2 .-= JacVert.B2_23[:,:,iz,ID] * rsC[2:3]
      #@views r3 .-= JacVert.B3_14[:,:,iz,ID] * rsC[[1,4]]
      for j = 1 : M2
        r2[j] -= JacVert.B2_23[j,1,iz,ID] * rs[sh + 2] + JacVert.B2_23[j,2,iz,ID] * rs[sh + 3] 
        r3[j] -= JacVert.B3_14[j,1,iz,ID] * rs[sh + 1] + JacVert.B3_14[j,2,iz,ID] * rs[sh + 4] 
      end  
      if iz > 1
        @views rsM = rs[sh - 1 : sh] 
        r1[1] -= JacVert.B1m_34[1,iz,ID] * rsM[1] + JacVert.B1m_34[2,iz,ID] * rsM[2]
      end
      if iz < nz
        @views rsP = rs[sh + 5 : sh + 6] 
        r1[M] -= JacVert.B1p_12[1,iz,ID] * rsP[1] + JacVert.B1p_12[2,iz,ID] * rsP[2]
      end
      @views ldivBlockAB(M,invfac,FacGrav,JacVert.A13[:,:,iz,ID],JacVert.A23[:,:,iz,ID],JacVert.A32[:,:,iz,ID],
        JacVert.SA[:,:,iz,ID],r1,r2,r3)
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
  @inbounds for iz = 1 : nz
    @inbounds for iv = 1 : 1
      @inbounds for k = 1 : M 
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
    @inbounds for iv = 2 : 3
      @inbounds for k = 2 : M - 1
        ii += 1
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end
    end
  end
  ivw = 3
  ivTh = 2
  @inbounds for iz = 1 : nz
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

#=
FillJacDGVert!(JacVert,U,DG,dz,fac,Phys,Param)
SchurBoundary!(JacVert)

b = ones(size(U))
@. b[:,:,:,5] *= 2
@. b[:,:,:,4] *= 3
ldivVertical!(JacVert,b)
=#

