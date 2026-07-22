@kernel inbounds = true function InterpolateRotKernel!(Rot,@Const(uC),@Const(vC),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(Inter),@Const(Glob))

  ID, iz   = @index(Local,  NTuple)
  _,Iz,IF = @index(Global,  NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform size(D,1)
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  RotLoc = @localmem eltype(Rot) (N,N,ColumnTilesDim)
  if Iz <= Nz
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    RotLoc[I,J,iz] = eltype(Rot)(0)
  end
  @synchronize

  if Iz <= Nz 
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    ind = Glob[ID,IF]
    tempx = ((dXdxI[1,1,1,ID,Iz,IF] + dXdxI[1,1,2,ID,Iz,IF]) * vC[Iz,ind] - 
      (dXdxI[1,2,1,ID,Iz,IF] + dXdxI[1,2,2,ID,Iz,IF]) * uC[Iz,ind])
    tempy = ((dXdxI[2,1,1,ID,Iz,IF] + dXdxI[2,1,2,ID,Iz,IF]) * vC[Iz,ind] - 
      (dXdxI[2,2,1,ID,Iz,IF] + dXdxI[2,2,2,ID,Iz,IF]) * uC[Iz,ind])
    for k = 1 : N
      @atomic :monotonic RotLoc[k,J,iz] += D[k,I] * tempx
      @atomic :monotonic RotLoc[I,k,iz] += D[k,J] * tempy
    end
  end
  @synchronize

  if Iz <= Nz 
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    temp = RotLoc[I,J,iz] / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    for j = 1 : size(Inter,1) 
      for i = 1 : size(Inter,2) 
        @atomic :monotonic Rot[i,j,Iz,IF] += Inter[i,j,1,I,J,1] * temp
      end    
    end    
  end    

end  

@kernel inbounds = true function InterpolateKernel!(cCell,@Const(c),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  I, J, K, iz   = @index(Local,  NTuple)
  _,_,_,Iz,IF = @index(Global,  NTuple)

  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]

  @uniform N = size(Inter,4)
  @uniform M = size(Inter,6)

  if Iz <= Nz
    cCell[I,J,K,Iz,IF] = eltype(cCell)(0)
    for kP = 1 : M
      iD = 0
      for jP = 1 : N
        for iP = 1 : N
          iD += 1  
          ind = Glob[iD,IF]
          cCell[I,J,K,Iz,IF] += Inter[I,J,K,iP,jP,kP] * c[kP,Iz,ind]
        end
      end
    end
  end
end

@kernel inbounds = true function InterpolateKernelNeu!(cCell,@Const(c),@Const(InterH),@Const(InterV),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  IP, KP, iz   = @index(Local,  NTuple)
  _,_,Iz,IF = @index(Global,  NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform DoFH = size(InterH,2)
  @uniform DoFV = size(InterV,2)

  if Iz <= Nz
    cCell[IP,KP,Iz,IF] = eltype(cCell)(0)  
    for iDoF = 1 : DoFH
      ind = Glob[iDoF,IF]  
      cLoc = eltype(cCell)(0)
      for k = 1 : DoFV  
        cLoc += InterV[KP,k] * c[k,Iz,ind] 
      end  
      cCell[IP,KP,Iz,IF] +=  InterH[IP,iDoF] * cLoc
    end
  end  
end

@kernel inbounds = true function InterpolateKernelCG!(cCell,@Const(c),@Const(InterH),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  IP, iz   = @index(Local,  NTuple)
  _,Iz,IF = @index(Global,  NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform DoFH = size(InterH,2)

  if Iz <= Nz
    cCell[IP,1,Iz,IF] = eltype(cCell)(0)
    for iDoF = 1 : DoFH
      ind = Glob[iDoF,IF]
      cCell[IP,1,Iz,IF] +=  InterH[IP,iDoF] * c[1,Iz,ind]
    end
  end
end

@kernel inbounds = true function InterpolateCGDim2Kernel!(cCell,@Const(c),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  I, J   = @index(Local,  NTuple)
  _,_,_,_,IF = @index(Global,  NTuple)


  ColumnTilesDim = @uniform @groupsize()[3]
  NF = @uniform @ndrange()[5]


  @uniform N = size(Inter,4)

  if IF <= NF
    cCell[I,J,IF] = eltype(cCell)(0)
    iD = 0
    for jP = 1 : N
      for iP = 1 : N
        iD += 1  
        ind = Glob[iD,IF]
        cCell[I,J,IF] += Inter[I,J,1,iP,jP,1] * c[1,1,ind]
      end
    end
  end
end

@kernel inbounds = true function InterpolateOrographyKernel!(cCell,@Const(c),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  I, J   = @index(Local,  NTuple)
  _,_,IF = @index(Global,  NTuple)


  ColumnTilesDim = @uniform @groupsize()[3]
  NF = @uniform @ndrange()[3]

  @uniform N = size(Inter,3)

  if IF <= NF
    cCell[I,J,IF] = eltype(cCell)(0)
    for jP = 1 : N
      for iP = 1 : N
         cCell[I,J,IF] += Inter[I,J,1,iP,jP,1] * c[I,J,IF]
      end
    end
  end
end

@kernel inbounds = true function InterpolateRhoKernel!(cCell,@Const(c),@Const(RhoC),@Const(InterH),
  @Const(InterV),@Const(Glob), ::Val{BANK}=Val(1)) where BANK

  IP, KP, iz   = @index(Local,  NTuple)
  _,_,Iz,IF = @index(Global,  NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform DoFH = size(InterH,2)
  @uniform DoFV = size(InterV,2)

  if Iz <= Nz
    cCellLoc = eltype(cCell)(0)
    for iDoF = 1 : DoFH
      ind = Glob[iDoF,IF]
      cLoc = eltype(cCell)(0)
      for k = 1 : DoFV
        cLoc += InterV[KP,k] * c[k,Iz,ind] / RhoC[k,Iz,ind]
      end
      cCellLoc +=  InterH[IP,iDoF] * cLoc
    end
    cCell[IP,KP,Iz,IF] = cCellLoc 
  end
end

@kernel inbounds = true function InterpolateWBKernel!(cCell,@Const(u),@Const(v),@Const(w),
 @Const(InterH),@Const(InterV),@Const(dXdxI),@Const(Glob))
  IP, KP, iz   = @index(Local,  NTuple)
  _,_,Iz,IF = @index(Global,  NTuple)

  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform DoFH = size(InterH,2)
  @uniform DoFV = size(InterV,2)

  if Iz == 1
    cCellLoc = eltype(cCell)(0)
    for iDoF = 1 : DoFH
      ind = Glob[iDoF,IF]
      cLoc = eltype(cCell)(0)
      for k = 1 : DoFV
        w0 = -(u[1,Iz,ind] * dXdxI[3,1,1,iDoF,Iz,IF] +
          v[1,Iz,ind] * dXdxI[3,2,1,iDoF,Iz,IF]) / dXdxI[3,3,1,iDoF,Iz,IF]  
        cLoc += eltype(cCell)(0.5) * InterV[KP,k] * (w[k,Iz,ind] + w0) 
      end
      cCellLoc +=  InterH[IP,iDoF] * cLoc
    end
    cCell[IP,KP,Iz,IF] = cCellLoc 
  elseif Iz <= Nz
    cCellLoc = eltype(cCell)(0)
    for iDoF = 1 : DoFH
      ind = Glob[iDoF,IF]
      cLoc = eltype(cCell)(0)
      for k = 1 : DoFV
        cLoc += eltype(cCell)(0.5) * InterV[KP,k] * (w[k,Iz,ind] + w[k,Iz-1,ind]) 
      end
      cCellLoc +=  InterH[IP,iDoF] * cLoc
    end
    cCell[IP,KP,Iz,IF] = cCellLoc 
  end
end

@kernel inbounds = true function InterpolateThEKernel!(cCell,@Const(RhoTh),@Const(Rho),@Const(RhoV),@Const(RhoC),@Const(Inter),@Const(Glob),@Const(Phys))
  I, J, K, iz   = @index(Local,  NTuple)
  _,_,_,Iz,IF = @index(Global,  NTuple)


  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]

  @uniform N = size(Inter,4)
  @uniform M = size(Inter,6)

  if Iz <= Nz
    cCell[I,J,K,Iz,IF] = eltype(cCell)(0)
    for kP = 1 : M
      iD = 0
      for jP = 1 : N
        for iP = 1 : N
          iD += 1
          ind = Glob[iD,IF]
          cLoc = Thermodynamics.fThE(Rho[kP,Iz,ind],RhoV[kP,Iz,ind],RhoC[kP,Iz,ind],
            RhoTh[Iz,kP,ind],Phys)
          cCell[I,J,K,Iz,IF] += Inter[I,J,K,iP,jP,kP] * cLoc
        end
      end
    end
  end    
end     

function InterpolateVortGPU!(Vort,U,Inter,FE::FiniteElements.CGElement,Metric,NumberThreadGPU,::Grids.Quad)

  backend = get_backend(U)
  FT = eltype(Vort)

  OrdPrint = size(Inter,1)
  DoF = FE.DoF
  Glob = FE.Glob
  NF = size(Glob,2)
  Nz = size(U,2)

# Ranges
  NzG = min(div(256,DoF),Nz)
  group = (DoF, NzG, 1)
  ndrange = (DoF, Nz, NF)

  @views uC = U[:,:,2]
  @views vC = U[:,:,3]
  @views dXdxI = Metric.dXdxI
  @views J = Metric.J
  D = FE.DS
  KInterpolateRotKernel! = InterpolateRotKernel!(backend,group)
  @. Vort = 0
  KInterpolateRotKernel!(Vort,uC,vC,D,dXdxI,J,Inter,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

function InterpolateVortGPU!(cCell,U,FE::FiniteElements.DGElement,Metric,NumberThreadGPU,::Grids.Quad)
  backend = get_backend(U)
  FT = eltype(U)
  Nz = size(U,2)
  NF = size(FE.Glob,2)
  NV = size(U,4)
  DoF = FE.DoF
  DoFE = FE.DoFE
  N = FE.OrdPoly + 1 
  M = FE.OrdPolyZ + 1 
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  KVorticityQuadKernel! = VorticityQuadKernel!(backend,group)
  Rot = KernelAbstractions.zeros(backend,FT,M,Nz,FE.NumI)
  KVorticityQuadKernel!(Rot,U,Metric.dXdxI,Metric.J,Metric.Rotate,FE.DS,FE.DSZ,FE.Glob,
    Val(N),Val(M) ;ndrange=ndrange)
  InterpolateGPU!(cCell,Rot,FE)
end

function InterpolateGPU!(cCell,c,FE::FiniteElements.DGElement)

  backend = get_backend(c)
  FT = eltype(c)

  InterH = FE.InterOutputH
  InterV = FE.InterOutputV
  Glob = FE.Glob 
  NumCellH = size(InterH,1)
  NumCellV = size(InterV,1)
  NF = size(Glob,2)
  Nz = size(c,2)
# Ranges
  NzG = min(div(256,NumCellH*NumCellV),Nz)
  group = (NumCellH,NumCellV,  NzG, 1)
  ndrange = (NumCellH,NumCellV, Nz, NF)
  KInterpolateKernel! = InterpolateKernelNeu!(backend,group)
  KInterpolateKernel!(cCell,c,InterH,InterV,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

function InterpolateGPU!(cCell,c,FE::FiniteElements.CGQuad)

  backend = get_backend(c)
  FT = eltype(c)

  InterH = FE.InterOutputH
  Glob = FE.Glob
  NumCellH = size(InterH,1)
  NF = size(Glob,2)
  Nz = size(c,2)
# Ranges
  NzG = min(div(256,NumCellH),Nz)
  group = (NumCellH,  NzG, 1)
  ndrange = (NumCellH, Nz, NF)
  KInterpolateKernel! = InterpolateKernelCG!(backend,group)
  KInterpolateKernel!(cCell,c,InterH,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end


function InterpolateCGDim2GPU!(cCell,c,Inter,Glob)

  backend = get_backend(c)
  FT = eltype(c)

  OrdPrint = size(Inter,1)
  OrdPrintZ = size(Inter,3)
  NF = size(cCell,3)
  Nz = size(c,2)
# Ranges
  NzG = min(div(256,OrdPrint*OrdPrint*OrdPrintZ),Nz)
  group = (OrdPrint, OrdPrint, OrdPrintZ,  NzG, 1)
  ndrange = (OrdPrint, OrdPrint, OrdPrintZ, Nz, NF)

  KInterpolateCGDim2Kernel! = InterpolateCGDim2Kernel!(backend,group)
  KInterpolateCGDim2Kernel!(cCell,c,Inter,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

function InterpolateRhoGPU!(cCell,c,Rho,FE)

  backend = get_backend(c)
  FT = eltype(c)

  InterH = FE.InterOutputH
  InterV = FE.InterOutputV
  Glob = FE.Glob
  NumCellH = size(InterH,1)
  NumCellV = size(InterV,1)
  NF = size(Glob,2)
  Nz = size(c,2)
# Ranges
  NzG = min(div(256,NumCellH*NumCellV),Nz)
  group = (NumCellH,NumCellV,  NzG, 1)
  ndrange = (NumCellH,NumCellV, Nz, NF)

  KInterpolateRhoKernel! = InterpolateRhoKernel!(backend,group)
  KInterpolateRhoKernel!(cCell,c,Rho,InterH,InterV,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

function InterpolateWBGPU!(cCell,u,v,w,dXdxI,FE)

  backend = get_backend(w)
  FT = eltype(w)

  InterH = FE.InterOutputH
  InterV = FE.InterOutputV
  Glob = FE.Glob
  NumCellH = size(InterH,1)
  NumCellV = size(InterV,1)
  NF = size(Glob,2)
  Nz = size(u,2)
# Ranges
  NzG = min(div(256,NumCellH*NumCellV),Nz)
  group = (NumCellH,NumCellV,  NzG, 1)
  ndrange = (NumCellH,NumCellV, Nz, NF)

  KInterpolateWBKernel! = InterpolateWBKernel!(backend,group)
  KInterpolateWBKernel!(cCell,u,v,w,InterH,InterV,dXdxI,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function InterpolateThEGPU!(cCell,RhoTh,Rho,RhoV,RhoC,Inter,Glob,Phys)

  backend = get_backend(cCell)
  FT = eltype(cCell)

  OrdPrint = size(Inter,1)
  OrdPrintZ = size(Inter,3)
  NF = size(Glob,2)
  Nz = size(RhoTh,2)
# Ranges
  NzG = min(div(256,OrdPrint*OrdPrint*OrdPrintZ),Nz)
  group = (OrdPrint, OrdPrint, OrdPrintZ,  NzG, 1)
  ndrange = (OrdPrint, OrdPrint, OrdPrintZ, Nz, NF)

  KInterpolateThEKernel! = InterpolateThEKernel!(backend,group)
  KInterpolateThEKernel!(cCell,RhoTh,Rho,RhoV,RhoC,Inter,Glob,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end  



@kernel inbounds = true function VorticityQuadKernel!(
    Rot, @Const(V), @Const(dXdxI), @Const(JJ), @Const(Rotate), @Const(D),
    @Const(DZ), @Const(Glob),
    ::Val{N}, ::Val{M}, ) where {N, M}

  I, J, K      = @index(Local, NTuple)
  _, _, _, Iz, IF  = @index(Global, NTuple)

  NZ       = @uniform @ndrange()[4]
  NF       = @uniform @ndrange()[5]

  VLoc     = @localmem eltype(Rot) (N, N, M, 3, 3)
  # 2 directions × 3 components
  RotLoc = @localmem eltype(Rot) (N, N, M)


  if Iz <= NZ && IF <= NF

    # ---- load phase ----
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    rho = V[K, Iz, ind, 1]
    vSp = @inbounds SVector{3}(V[K,Iz,ind,2] / rho, V[K,Iz,ind,3] / rho, V[K,Iz,ind,4] / rho)

    R = SMatrix{3, 3}(
        Rotate[1, 1, ID, IF], Rotate[1, 2, ID, IF], Rotate[1, 3, ID, IF],
        Rotate[2, 1, ID, IF], Rotate[2, 2, ID, IF], Rotate[2, 3, ID, IF],
        Rotate[3, 1, ID, IF], Rotate[3, 2, ID, IF], Rotate[3, 3, ID, IF]
      )

    vCa = R * vSp

    m  = SVector{3}(dXdxI[1,1,K,ID,Iz,IF], dXdxI[1,2,K,ID,Iz,IF], dXdxI[1,3,K,ID,Iz,IF])
    u = cross(m,vCa)
    VLoc[I,J,K,1,1] = u[1]
    VLoc[I,J,K,2,1] = u[2]
    VLoc[I,J,K,3,1] = u[3]
    m  = SVector{3}(dXdxI[2,1,K,ID,Iz,IF], dXdxI[2,2,K,ID,Iz,IF], dXdxI[2,3,K,ID,Iz,IF])
    u = cross(m,vCa)
    VLoc[I,J,K,1,2] = u[1]
    VLoc[I,J,K,2,2] = u[2]
    VLoc[I,J,K,3,2] = u[3]
    m  = SVector{3}(dXdxI[3,1,K,ID,Iz,IF], dXdxI[3,2,K,ID,Iz,IF], dXdxI[3,3,K,ID,Iz,IF])
    u = cross(m,vCa)
    VLoc[I,J,K,1,3] = u[1]
    VLoc[I,J,K,2,3] = u[2]
    VLoc[I,J,K,3,3] = u[3]
  end  

  @synchronize

  if Iz <= NZ && IF <= NF

    ID = I + (J - 1) * N
    ind = Glob[ID, IF]

    dU = D[I,1] * VLoc[1,J,K,1,1] + D[J,1] * VLoc[I,1,K,1,2]
    dV = D[I,1] * VLoc[1,J,K,2,1] + D[J,1] * VLoc[I,1,K,2,2] 
    dW = D[I,1] * VLoc[1,J,K,3,1] + D[J,1] * VLoc[I,1,K,3,2]
    @unroll for l = 2 : N
      dU += D[I,l] * VLoc[l,J,K,1,1] + D[J,l] * VLoc[I,l,K,1,2]
      dV += D[I,l] * VLoc[l,J,K,2,1] + D[J,l] * VLoc[I,l,K,2,2]
      dW += D[I,l] * VLoc[l,J,K,3,1] + D[J,l] * VLoc[I,l,K,3,2]
    end  
    @unroll for l = 1 : M
      dU += DZ[K,l] * VLoc[I,J,l,1,3]
      dV += DZ[K,l] * VLoc[I,J,l,2,3]
      dW += DZ[K,l] * VLoc[I,J,l,3,3]
    end  
    vCa = @inbounds SVector{3}(dU,dV,dW)

    R = SMatrix{3, 3}(
        Rotate[1, 1, ID, IF], Rotate[2, 1, ID, IF], Rotate[3, 1, ID, IF],
        Rotate[1, 2, ID, IF], Rotate[2, 2, ID, IF], Rotate[3, 2, ID, IF],
        Rotate[1, 3, ID, IF], Rotate[2, 3, ID, IF], Rotate[3, 3, ID, IF]
      )

    vSp = R * vCa
    JLoc = JJ[ID,K,Iz,IF]
    Rot[K, Iz, ind] = vSp[3] * JLoc
end  

@kernel inbounds = true function InterpolateKernel!(cCell,@Const(c),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  I, J, K, iz   = @index(Local,  NTuple)
  _,_,_,Iz,IF = @index(Global,  NTuple)

  end  
end  
