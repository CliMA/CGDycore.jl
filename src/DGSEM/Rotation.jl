# Rotation of the momentum from Cartesian components to the local frame
# (zonal, meridional, radial) of each column, used for the output
function RotateCart2Sph!(V,DG,Metric,NumberThreadGPU)
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

function ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  backend = get_backend(F)
  M = DG.V.OrdPolyZ + 1
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



