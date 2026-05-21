function TendVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
  backend = get_backend(V)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  NFG = min(div(NumberThreadGPU,DoF),NF)
  group = (DoF,NFG)
  ndrange = (DoF,NF)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  KVSp2VCart3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
end  

function StateVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end

function TendVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
  backend = get_backend(V)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  NFG = min(div(NumberThreadGPU,DoF),NF)
  group = (DoF,NFG)
  ndrange = (DoF,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
end  

function StateVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end  
function StateVCart2VSpScale!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end  

function StateVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  DoFG = min(div(NumberThreadGPU,M*Nz),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  KVSp2VCart3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
end  

function TendVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
end

function StateVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  DoFG = min(div(NumberThreadGPU,M*Nz),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
end  

function StateVCart2VSpScale!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  DoFG = min(div(NumberThreadGPU,M*Nz),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KVCart2VSp3ScaleKernel! = VCart2VSp3ScaleKernel!(backend,group)
  KVCart2VSp3ScaleKernel!(V,Metric.Rotate,Metric.J,DG.Glob;ndrange=ndrange)
end

function TendVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
end  

@kernel inbounds = true function VSp2VCartKernel!(VCart,@Const(VSp),@Const(Rotate))

  K,iQ,  = @index(Local, NTuple) 
  _,Iz,IF = @index(Global, NTuple)
  
  M = @uniform @ndrange()[1]
  DoF = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if IF <= NF
    iD = iQ + (IF - 1) * DoF
    VCart[1,1,iD,1] = VSp[1,1,iD,1]
    VCart[1,1,iD,2] = Rotate[1,1,1,iQ,1,IF] * VSp[1,1,iD,2] +
      Rotate[2,1,1,iQ,1,IF] * VSp[1,1,iD,3]
    VCart[1,1,iD,3] = Rotate[1,2,1,iQ,1,IF] * VSp[1,1,iD,2] +
      Rotate[2,2,1,iQ,1,IF] * VSp[1,1,iD,3]
    VCart[1,1,iD,4] = Rotate[1,3,1,iQ,1,IF] * VSp[1,1,iD,2] +
      Rotate[2,3,1,iQ,1,IF] * VSp[1,1,iD,3]
  end
end

@kernel inbounds = true function VCart2VSpKernel!(VSp,@Const(VCart),@Const(Rotate),@Const(J))

  iQ,  = @index(Local, NTuple)
  _,IF = @index(Global, NTuple)

  DoF = @uniform @ndrange()[1]
  NF = @uniform @ndrange()[2]

  if IF <= NF
    iD = iQ + (IF - 1) * DoF
    VSp[1,1,iD,1] = VCart[1,1,iD,1] / J[iQ,1,1,IF]
    VSp[1,1,iD,2] = (Rotate[1,1,1,iQ,1,IF] * VCart[1,1,iD,2] +
      Rotate[1,2,1,iQ,1,IF] * VCart[1,1,iD,3] +
      Rotate[1,3,1,iQ,1,IF] * VCart[1,1,iD,4]) / J[iQ,1,1,IF]
    VSp[1,1,iD,3] = (Rotate[2,1,1,iQ,1,IF] * VCart[1,1,iD,2] +
      Rotate[2,2,1,iQ,1,IF] * VCart[1,1,iD,3] +
      Rotate[2,3,1,iQ,1,IF] * VCart[1,1,iD,4]) / J[iQ,1,1,IF]
  end
end

@kernel inbounds = true function VSp2VCart3Kernel!(V, @Const(Rotate), @Const(Glob))

  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    # Load the 3x3 slice directly into a static matrix (stored in registers)
    R = SMatrix{3, 3}(
      Rotate[1, 1, ID, IF], Rotate[1, 2, ID, IF], Rotate[1, 3, ID, IF],
      Rotate[2, 1, ID, IF], Rotate[2, 2, ID, IF], Rotate[2, 3, ID, IF],
      Rotate[3, 1, ID, IF], Rotate[3, 2, ID, IF], Rotate[3, 3, ID, IF]
    )

    ind = Glob[ID, IF]

    # Load velocity component into a static vector
    v = @inbounds SVector{3}(V[K, Iz, ind, 2], V[K, Iz, ind, 3], V[K, Iz, ind, 4])

    # Matrix-vector multiplication (completely optimized in registers)
    v_rot = R * v

    # Write back to global memory
    V[K, Iz, ind, 2] = v_rot[1]
    V[K, Iz, ind, 3] = v_rot[2]
    V[K, Iz, ind, 4] = v_rot[3]
  end
end

@kernel inbounds = true function VCart2VSp3Kernel!(V, @Const(Rotate), @Const(Glob))

  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    # Load the 3x3 slice directly into a static matrix (stored in registers)
    R = SMatrix{3, 3}(
      Rotate[1, 1, ID, IF], Rotate[2, 1, ID, IF], Rotate[3, 1, ID, IF],
      Rotate[1, 2, ID, IF], Rotate[2, 2, ID, IF], Rotate[3, 2, ID, IF],
      Rotate[1, 3, ID, IF], Rotate[2, 3, ID, IF], Rotate[3, 3, ID, IF]
    )

    ind = Glob[ID, IF]

    # Load velocity component into a static vector
    v = SVector{3}(V[K, Iz, ind, 2], V[K, Iz, ind, 3], V[K, Iz, ind, 4])

    # Matrix-vector multiplication (completely optimized in registers)
    v_rot = R * v

    # Write back to global memory
    V[K, Iz, ind, 2] = v_rot[1]
    V[K, Iz, ind, 3] = v_rot[2]
    V[K, Iz, ind, 4] = v_rot[3]
  end
end

@kernel inbounds = true function VCart2VSp3ScaleKernel!(V,@Const(Rotate),@Const(J),@Const(Glob))

  K,Iz,ID,IF = @index(Global, NTuple)

  ND = @uniform @ndrange()[3]

  if ID <= ND
    # Load the 3x3 slice directly into a static matrix (stored in registers)
    R = SMatrix{3, 3}(
      Rotate[1, 1, ID, IF], Rotate[2, 1, ID, IF], Rotate[3, 1, ID, IF],
      Rotate[1, 2, ID, IF], Rotate[2, 2, ID, IF], Rotate[3, 2, ID, IF],
      Rotate[1, 3, ID, IF], Rotate[2, 3, ID, IF], Rotate[3, 3, ID, IF]
    )

    ind = Glob[ID, IF]

    # Load velocity component into a static vector
    v = SVector{3}(V[K, Iz, ind, 2], V[K, Iz, ind, 3], V[K, Iz, ind, 4])

    # Matrix-vector multiplication (completely optimized in registers)
    v_rot = R * v

    # Write back to global memory
    JLoc = J[K,ID,Iz,IF]  
    V[K,Iz,ind,1] *= JLoc
    V[K,Iz,ind,2] = v_rot[1] * JLoc
    V[K,Iz,ind,3] = v_rot[2] * JLoc
    V[K,Iz,ind,4] = v_rot[3] * JLoc
    V[K,Iz,ind,5] *= JLoc
  end
end

function ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  backend = get_backend(F)
  DoF = DG.DoF
  M = DG.OrdPolyZ + 1
  Nz = Grid.nz
  NF = Grid.NumFaces
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (M,Nz,DoFG,1)
  ndrange = (M,Nz,DoF,NF)
  KScaleMassMatrixKernel! = ScaleMassMatrixKernel!(backend,group)
  KScaleMassMatrixKernel!(F,Metric.J,DG.Glob,Val(NV);ndrange=ndrange)
end  

@kernel inbounds = true function ScaleMassMatrixKernel!(F,@Const(J),@Const(Glob), 
  ::Val{NUMV}) where {NUMV}

  _,_,iD,  = @index(Local, NTuple)
  K,Iz,ID,IF = @index(Global, NTuple)

  DoF = @uniform @ndrange()[3]

  if ID <= DoF
    ind = Glob[ID,IF]
    JLoc = J[K,ID,Iz,IF]
    @unroll for iv = 1 : NUMV
      F[K,Iz,ind,iv] *= JLoc
    end  
  end
end



