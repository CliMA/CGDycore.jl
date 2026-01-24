function TendVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KVSp2VCart3Kernel! = VSp2VCart3Kernel!(backend,group)
  KVSp2VCart3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
end  

function StateVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end

function TendVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
end  

function StateVCart2VSp!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityC)
end  

function StateVSp2VCart!(V,DG,Metric,NumberThreadGPU,::Examples.VelocityS)
  backend = get_backend(V)
  M = size(V,1)
  Nz = size(V,2)
  NF = size(DG.Glob,2)
  DoF = DG.DoF
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
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
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KVCart2VSp3Kernel! = VCart2VSp3Kernel!(backend,group)
  KVCart2VSp3Kernel!(V,Metric.Rotate,DG.Glob;ndrange=ndrange)
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

@kernel inbounds = true function VSp2VCart3Kernel!(V,@Const(Rotate),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  DoF = @uniform @ndrange()[3]

  if ID <= DoF
    ind = Glob[ID,IF]
    v1 = V[K,Iz,ind,1]
    v2 = V[K,Iz,ind,2]
    v3 = V[K,Iz,ind,3]
    V[K,Iz,ind,1] = Rotate[1,1,K,ID,Iz,IF] * v1 +
      Rotate[2,1,K,ID,Iz,IF] * v2 +
      Rotate[3,1,K,ID,Iz,IF] * v3
    V[K,Iz,ind,2] = Rotate[1,2,1,ID,Iz,IF] * v1 +
      Rotate[2,2,K,ID,Iz,IF] * v2 +
      Rotate[3,2,K,ID,Iz,IF] * v3
    V[K,Iz,ind,3] = Rotate[1,3,K,ID,Iz,IF] * v1 +
      Rotate[2,3,K,ID,Iz,IF] * v2 +
      Rotate[3,3,K,ID,Iz,IF] * v3
  end  
end


@kernel inbounds = true function VCart2VSp3Kernel!(V,@Const(Rotate),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  DoF = @uniform @ndrange()[3]
  M = @uniform @ndrange()[2]
  Nz = @uniform @ndrange()[1]

  if ID <= DoF
    ind = Glob[ID,IF]
    v1 = V[K,Iz,ind,1]
    v2 = V[K,Iz,ind,2]
    v3 = V[K,Iz,ind,3]
    V[K,Iz,ind,1] = Rotate[1,1,K,ID,Iz,IF] * v1 +
      Rotate[1,2,K,ID,Iz,IF] * v2 +
      Rotate[1,3,K,ID,Iz,IF] * v3
    V[K,Iz,ind,2] = Rotate[2,1,K,ID,Iz,IF] * v1 +
      Rotate[2,2,K,ID,Iz,IF] * v2 +
      Rotate[2,3,K,ID,Iz,IF] * v3
#   if Nz == 1 && M == 2
#     V[K,Iz,ind,3] = eltype(V)(0)  
#   else  
      V[K,Iz,ind,3] = Rotate[3,1,K,ID,Iz,IF] * v1 +
        Rotate[3,2,K,ID,Iz,IF] * v2 +
        Rotate[3,3,K,ID,Iz,IF] * v3
#   end  
#   if Iz == Nz && K == M
#     V[K,Iz,ind,3] = eltype(V)(0)  
#   end  
  end  
end

function ScaleMassMatrix!(F,DG,Metric,Grid,NumberThreadGPU,NV)

  backend = get_backend(F)
  DoF = DG.DoF
  M = DG.OrdPolyZ + 1
  Nz = Grid.nz
  NF = Grid.NumFaces
  DoFG = min(div(NumberThreadGPU,Nz*M),DoF)
  group = (Nz,M,DoFG,1)
  ndrange = (Nz,M,DoF,NF)
  KScaleMassMatrixKernel! = ScaleMassMatrixKernel!(backend,group)
  KScaleMassMatrixKernel!(F,Metric.J,DG.Glob,Val(NV);ndrange=ndrange)
end  

@kernel inbounds = true function ScaleMassMatrixKernel!(F,@Const(J),@Const(Glob), 
  ::Val{NUMV}) where {NUMV}

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  DoF = @uniform @ndrange()[3]

  if ID <= DoF
    ind = Glob[ID,IF]
    @unroll for iv = 1 : NUMV
      F[K,Iz,ind,iv] /= J[ID,K,Iz,IF]
    end  
  end
end



