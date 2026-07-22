function Vorticity(Vort,U,DG::FiniteElements.DGElement,Metric,NumberThreadGPU,::Grids.Quad)
  backend = get_backend(Vort)
  Nz = size(U,2)
  NF = size(DG.Glob,2)
  NV = size(U,4)
  DoF = DG.DoF
  DoFE = DG.DoFE
  N = DG.OrdPoly + 1 
  M = DG.OrdPolyZ + 1 
  group = (N,N,M,1,1)
  ndrange = (N,N,M,Nz,NF)
  KVorticityQuadKernel! = VorticityQuadKernel!(backend,group)
  KVorticityQuadKernel!(Vort,U,Metric.dXdxI,Metric.J,Metric.Rotate,DG.DS,DG.DSZ,DG.Glob,
    Val(N),Val(M) ;ndrange=ndrange)
end


@kernel inbounds = true function VorticityQuadKernel!(
    Vort, @Const(V), @Const(dXdxI), @Const(JJ), @Const(Rotate), @Const(D),
    @Const(DZ), @Const(Glob),
    ::Val{N}, ::Val{M}, ) where {N, M}

  I, J, K      = @index(Local, NTuple)
  _, _, _, Iz, IF  = @index(Global, NTuple)

  NZ       = @uniform @ndrange()[4]
  NF       = @uniform @ndrange()[5]

  VLoc     = @localmem eltype(Vort) (N, N, M, 3, 3)
  # 2 directions × 3 components
  dXdxILoc = @localmem eltype(Vort) (3, 3, N, N, M)
  DLoc = @private eltype(Vort) (N, N)
  DZLoc = @private eltype(Vort) (M, M)

  DLoc .= D
  DZLoc .= DZ

  if Iz <= NZ && IF <= NF

    # ---- load phase ----
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]
    rho = V[K, Iz, ind, 2]
    vSp = @inbounds SVector{3}(V[K, Iz, ind, 2] / rho, V[K, Iz, ind, 3] / rho, V[K, Iz, ind, 4] / rho)

    R = SMatrix{3, 3}(
        Rotate[1, 1, ID, IF], Rotate[1, 2, ID, IF], Rotate[1, 3, ID, IF],
        Rotate[2, 1, ID, IF], Rotate[2, 2, ID, IF], Rotate[2, 3, ID, IF],
        Rotate[3, 1, ID, IF], Rotate[3, 2, ID, IF], Rotate[3, 3, ID, IF]
      )

    vCa = R * vSp

    m  = SVector{3}(dXdxI[1, 1, ID, Iz, IF], dXdxI[1, 2, ID, Iz, IF], dXdxI[1, 3, ID, Iz, IF])
    u = cross(m,vCa)
    VLoc[I,J,K,1,1] = u[1]
    VLoc[I,J,K,2,1] = u[2]
    VLoc[I,J,K,3,1] = u[3]
    m  = SVector{3}(dXdxI[2, 1, ID, Iz, IF], dXdxI[2, 2, ID, Iz, IF], dXdxI[2, 3, ID, Iz, IF])
    u = cross(m,vCa)
    VLoc[I,J,K,1,2] = u[1]
    VLoc[I,J,K,2,2] = u[2]
    VLoc[I,J,K,3,2] = u[3]
    m  = SVector{3}(dXdxI[3, 1, ID, Iz, IF], dXdxI[3, 2, ID, Iz, IF], dXdxI[3, 3, ID, Iz, IF])
    u = cross(m,vCa)
    VLoc[I,J,K,1,3] = u[1]
    VLoc[I,J,K,2,3] = u[2]
    VLoc[I,J,K,3,3] = u[3]
  end  

  @synchronize

  if Iz <= NZ && IF <= NF

    # ---- load phase ----
    ID = I + (J - 1) * N
    ind = Glob[ID, IF]

    dU = DLoc[I,1] * VLoc[1,J,K,1,1] + DLoc[J,1] * VLoc[I,1,K,1,2]
    dV = DLoc[I,1] * VLoc[1,J,K,2,1] + DLoc[J,1] * VLoc[I,1,K,2,2] 
    dW = DLoc[I,1] * VLoc[1,J,K,3,1] + DLoc[J,1] * VLoc[I,1,K,3,2]
    @unroll for l = 2 : N
      dU += DLoc[I,l] * VLoc[l,J,K,1,1] + DLoc[J,l] * VLoc[I,l,K,1,2]
      dV += DLoc[I,l] * VLoc[l,J,K,2,1] + DLoc[J,l] * VLoc[I,l,K,2,2]
      dW += DLoc[I,l] * VLoc[l,J,K,3,1] + DLoc[J,l] * VLoc[I,l,K,3,2]
    end  
    @unroll for l = 1 : M
      dU += DZLoc[I,l] * VLoc[I,J,l,1,3]
      dV += DZLoc[I,l] * VLoc[I,J,l,2,3]
      dW += DZLoc[I,l] * VLoc[I,J,l,3,3]
    end  
    vCa = @inbounds SVector{3}(dU,dV,dW)

    R = SMatrix{3, 3}(
        Rotate[1, 1, ID, IF], Rotate[2, 1, ID, IF], Rotate[3, 1, ID, IF],
        Rotate[1, 2, ID, IF], Rotate[2, 2, ID, IF], Rotate[3, 2, ID, IF],
        Rotate[1, 3, ID, IF], Rotate[2, 3, ID, IF], Rotate[3, 3, ID, IF]
      )

    vSp = R * vCa
    JLoc = JJ[ID,1,Iz,IF]
    Vort[K, Iz, ind, 1] = vSp[1] / JLoc
    Vort[K, Iz, ind, 2] = vSp[2] / JLoc
    Vort[K, Iz, ind, 3] = vSp[3] / JLoc
  end  
end  
