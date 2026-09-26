function TendVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  ND = DG.NumI
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG)
  ndrange = (M,Nz,ND)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  KVSp2VCart3Kernel!(V,Metric.Rotate;ndrange=ndrange)
end  

function StateVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end

function TendVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  ND = DG.NumI
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG)
  ndrange = (M,Nz,ND)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(V,Metric.Rotate;ndrange=ndrange)
end  

function StateVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end  
function StateVCart2VSpScale!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end  

function StateVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  ND = DG.NumI
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG)
  ndrange = (M,Nz,ND)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  KVSp2VCart3Kernel!(V,Metric.Rotate;ndrange=ndrange)
end  

function TendVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
end

function StateVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  ND = DG.NumI
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG)
  ndrange = (M,Nz,ND)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(V,Metric.Rotate;ndrange=ndrange)
end  

function StateVCart2VSpScale!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  ND = DG.NumI
  NDG = min(div(NumberThreadGPU,M*Nz),ND)
  group = (M,Nz,NDG)
  ndrange = (M,Nz,ND)
  KVCart2VSp3ScaleKernel! = VCart2VSp3ScaleKernel!(backend,group)
  KVCart2VSp3ScaleKernel!(V,Metric.Rotate,Metric.J;ndrange=ndrange)
end

function TendVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
end  

@kernel inbounds = true function VSp2VCart3Kernel!(V, @Const(Rotate))

  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND

    v1 = V[K, Iz, ID, 2]
    v2 = V[K, Iz, ID, 3]
    v3 = V[K, Iz, ID, 4]

    V[K, Iz, ID, 2] = Rotate[1, 1, ID] * v1 + Rotate[2, 1, ID] * v2 + Rotate[3, 1, ID] * v3
    V[K, Iz, ID, 3] = Rotate[1, 2, ID] * v1 + Rotate[2, 2, ID] * v2 + Rotate[3, 2, ID] * v3
    V[K, Iz, ID, 4] = Rotate[1, 3, ID] * v1 + Rotate[2, 3, ID] * v2 + Rotate[3, 3, ID] * v3
  end
end



@kernel inbounds = true function VCart2VSp3Kernel!(V, @Const(Rotate))

  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND

    v1 = V[K, Iz, ID, 2]
    v2 = V[K, Iz, ID, 3]
    v3 = V[K, Iz, ID, 4]

    V[K, Iz, ID, 2] = Rotate[1, 1, ID] * v1 + Rotate[1, 2, ID] * v2 + Rotate[1, 3, ID] * v3
    V[K, Iz, ID, 3] = Rotate[2, 1, ID] * v1 + Rotate[2, 2, ID] * v2 + Rotate[2, 3, ID] * v3
    V[K, Iz, ID, 4] = Rotate[3, 1, ID] * v1 + Rotate[3, 2, ID] * v2 + Rotate[3, 3, ID] * v3
  end
end

@kernel inbounds = true function VCart2VSp3ScaleKernel!(V,@Const(Rotate),@Const(invJ))

  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    # Load the 3x3 slice directly into a static matrix (stored in registers)
    R = SMatrix{3, 3}(
      Rotate[1, 1, ID], Rotate[2, 1, ID], Rotate[3, 1, ID],
      Rotate[1, 2, ID], Rotate[2, 2, ID], Rotate[3, 2, ID],
      Rotate[1, 3, ID], Rotate[2, 3, ID], Rotate[3, 3, ID]
    )


    # Load velocity component into a static vector
    v = SVector{3}(V[K, Iz, ID, 2], V[K, Iz, ID, 3], V[K, Iz, ID, 4])

    # Matrix-vector multiplication (completely optimized in registers)
    v_rot = R * v

    # Write back to global memory
    invJLoc = invJ[K,Iz,ID]  
    V[K,Iz,ID,1] *= invJLoc
    V[K,Iz,ID,2] = v_rot[1] * invJLoc
    V[K,Iz,ID,3] = v_rot[2] * invJLoc
    V[K,Iz,ID,4] = v_rot[3] * invJLoc
    V[K,Iz,ID,5] *= invJLoc
  end
end

function ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  backend = get_backend(F)
  M = DG.OrdPolyZ + 1
  Nz = Grid.nz
  ND = DG.NumI
  NDG = min(div(NumberThreadGPU,Nz*M),ND)
  group = (M,Nz,NDG)
  ndrange = (M,Nz,ND)
  KScaleMassMatrixKernel! = ScaleMassMatrixKernel!(backend,group)
  KScaleMassMatrixKernel!(F,Metric.invJ,Val(NV);ndrange=ndrange)
end  

@kernel inbounds = true function ScaleMassMatrixKernel!(F,@Const(invJ), 
  ::Val{NUMV}) where {NUMV}

  K,Iz,ID = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    invJLoc = invJ[K,Iz,ID]
    @unroll for iv = 1 : NUMV
      F[K,Iz,ID,iv] *= invJLoc
    end  
  end
end



