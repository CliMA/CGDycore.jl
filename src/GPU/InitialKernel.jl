@kernel function uvwFunCKernel!(Profile,u,v,w,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]


  if Iz <= Nz
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,uP,vP,_ = Profile(xS,time)
    @inbounds u[Iz,ind] = uP
    @inbounds v[Iz,ind] = vP
  end
  if Iz <= Nz - 1
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,2,1,Iz,IF] + X[I,1,1,Iz+1,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,2,2,Iz,IF] + X[I,1,2,Iz+1,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,2,3,Iz,IF] + X[I,1,3,Iz+1,IF])
    xS = SVector{3}(x1, x2 ,x3)
    @inbounds _,_,_,w[Iz,ind] = Profile(xS,time)
  end
end

@kernel function RhouvFunCKernel!(Profile,Rho,u,v,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    @inbounds ind = Glob[ID,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,uP,vP,_ = Profile(xS,time)
    @inbounds Rho[Iz,ind] = RhoP
    @inbounds u[Iz,ind] = uP
    @inbounds v[Iz,ind] = vP
  end
end

@kernel function RhoFunCKernel!(Profile,Rho,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ = Profile(xS,time)
    @inbounds Rho[Iz,ind] = RhoP
  end
end

@kernel function TrFunCKernel!(Profile,Tr,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,TrP = Profile(xS,time)
    @inbounds Tr[Iz,ind] = RhoP * TrP
  end
end

@kernel function RhoThFunCKernel!(Profile,RhoTh,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,ThP = Profile(xS,time)
    @inbounds RhoTh[Iz,ind] = RhoP * ThP
  end
end

@kernel function RhoVFunCKernel!(Profile,RhoV,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_,_,QvP = Profile(xS,time)
    @inbounds RhoV[Iz,ind] = RhoP * QvP
  end
end

@kernel function RhoCFunCKernel!(Profile,RhoC,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    @inbounds ind = Glob[I,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_,_,_,QcP = Profile(xS,time)
    @inbounds RhoC[Iz,ind] = RhoP * QcP
  end
end

@kernel function ComputeFunFKernel!(Profile,w,time,@Const(Glob),@Const(X),Param)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz - 1 
    @inbounds ind = Glob[ID,IF]
    @inbounds x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    @inbounds x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    @inbounds x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,_,_,wP = Profile(xS,time)
    @inbounds w[Iz,ind] = wP
  end
end
