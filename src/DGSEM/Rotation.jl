@kernel inbounds = true function VSp2VCartKernel!(VCart,@Const(VSp),@Const(Rotate))

  K,iQ,  = @index(Local, NTuple) 
  _,Iz,IF = @index(Global, NTuple)
  
  M = @uniform @ndrange()[1]
  NQ = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if IF <= NF
    iD = iQ + (IF - 1) * NQ
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

  NQ = @uniform @ndrange()[1]
  NF = @uniform @ndrange()[2]

  if IF <= NF
    iD = iQ + (IF - 1) * NQ
    VSp[1,1,iD,1] = VCart[1,1,iD,1] / J[iQ,1,1,IF]
    VSp[1,1,iD,2] = (Rotate[1,1,1,iQ,1,IF] * VCart[1,1,iD,2] +
      Rotate[1,2,1,iQ,1,IF] * VCart[1,1,iD,3] +
      Rotate[1,3,1,iQ,1,IF] * VCart[1,1,iD,4]) / J[iQ,1,1,IF]
    VSp[1,1,iD,3] = (Rotate[2,1,1,iQ,1,IF] * VCart[1,1,iD,2] +
      Rotate[2,2,1,iQ,1,IF] * VCart[1,1,iD,3] +
      Rotate[2,3,1,iQ,1,IF] * VCart[1,1,iD,4]) / J[iQ,1,1,IF]
  end
end

@kernel inbounds = true function VSp2VCart3Kernel!(VCart,@Const(VSp),@Const(Rotate),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[3]

  if ID <= NQ
    ind = Glob[ID,IF]
    VCart[Iz,K,ind,1] = VSp[Iz,K,ind,1]
    VCart[Iz,K,ind,2] = Rotate[1,1,K,ID,Iz,IF] * VSp[Iz,K,ind,2] +
      Rotate[2,1,K,ID,Iz,IF] * VSp[Iz,K,ind,3] +
      Rotate[3,1,K,ID,Iz,IF] * VSp[Iz,K,ind,4]
    VCart[Iz,K,ind,3] = Rotate[1,2,1,ID,Iz,IF] * VSp[Iz,K,ind,2] +
      Rotate[2,2,K,ID,Iz,IF] * VSp[Iz,K,ind,3] +
      Rotate[3,2,K,ID,Iz,IF] * VSp[Iz,K,ind,4] 
    VCart[Iz,K,ind,4] = Rotate[1,3,K,ID,Iz,IF] * VSp[Iz,K,ind,2] +
      Rotate[2,3,K,ID,Iz,IF] * VSp[Iz,K,ind,3] +
      Rotate[3,3,K,ID,Iz,IF] * VSp[Iz,K,ind,4]
    VCart[Iz,K,ind,5] = VSp[Iz,K,ind,5]
  end  
end


@kernel inbounds = true function VCart2VSp3Kernel!(VSp,@Const(VCart),@Const(Rotate),@Const(Glob))

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[3]
  M = @uniform @ndrange()[2]
  Nz = @uniform @ndrange()[1]

  if ID <= NQ
    ind = Glob[ID,IF]
    VSp[Iz,K,ind,1] = VCart[Iz,K,ind,1] 
    VSp[Iz,K,ind,2] = (Rotate[1,1,K,ID,Iz,IF] * VCart[Iz,K,ind,2] +
      Rotate[1,2,K,ID,Iz,IF] * VCart[Iz,K,ind,3] +
      Rotate[1,3,K,ID,Iz,IF] * VCart[Iz,K,ind,4])
    VSp[Iz,K,ind,3] = (Rotate[2,1,K,ID,Iz,IF] * VCart[Iz,K,ind,2] +
      Rotate[2,2,K,ID,Iz,IF] * VCart[Iz,K,ind,3] +
      Rotate[2,3,K,ID,Iz,IF] * VCart[Iz,K,ind,4])
    if Nz == 1 && M == 2
      VSp[Iz,K,ind,4] = eltype(VSp)(0)  
    else  
      VSp[Iz,K,ind,4] = (Rotate[3,1,K,ID,Iz,IF] * VCart[Iz,K,ind,2] +
        Rotate[3,2,K,ID,Iz,IF] * VCart[Iz,K,ind,3] +
        Rotate[3,3,K,ID,Iz,IF] * VCart[Iz,K,ind,4])
    end  
    VSp[Iz,K,ind,5] = VCart[Iz,K,ind,5]
  end  
end

@kernel inbounds = true function ScaleMassMatrixKernel!(F,@Const(J),@Const(Glob), 
  ::Val{NUMV}) where {NUMV}

  _,_,iD,  = @index(Local, NTuple)
  Iz,K,ID,IF = @index(Global, NTuple)

  NQ = @uniform @ndrange()[3]

  if ID <= NQ
    ind = Glob[ID,IF]
    @unroll for iv = 1 : NUMV
      F[Iz,K,ind,iv] /= J[ID,K,Iz,IF]
    end  
  end
end



