@kernel inbounds = true function uvwFunCKernel!(Profile,u,v,w,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]


  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,uP,vP,_ = Profile(xS,time)
    u[1,Iz,ind] = uP
    v[1,Iz,ind] = vP
  end
  if Iz <= Nz - 1
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,2,1,Iz,IF] + X[I,1,1,Iz+1,IF])
    x2 = eltype(X)(0.5) * (X[I,2,2,Iz,IF] + X[I,1,2,Iz+1,IF])
    x3 = eltype(X)(0.5) * (X[I,2,3,Iz,IF] + X[I,1,3,Iz+1,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,_,_,w[1,Iz,ind] = Profile(xS,time)
  end
end

@kernel inbounds = true function uvwFunCDGKernel!(Profile,u,v,w,time,@Const(Glob),@Const(X),Param,Phys)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  M = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]


  if Iz <= Nz
    ind = Glob[ID,IF]
    x1 = X[ID,K,1,Iz,IF] 
    x2 = X[ID,K,2,Iz,IF]
    x3 = X[ID,K,3,Iz,IF] 
    xS = SVector{3}(x1, x2 ,x3)
    _,uP,vP,wP = Profile(xS,time)
    u[K,Iz,ind] = uP
    v[K,Iz,ind] = vP
    w[K,Iz,ind] = wP
  end
end

@kernel inbounds = true function uvFunCDGKernel!(Profile,u,v,time,@Const(Glob),@Const(X),Param,Phys)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  M = @uniform @groupsize()[1]
  N = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]


  if Iz <= Nz
    ind = Glob[ID,IF]
    x1 = X[ID,K,1,Iz,IF] 
    x2 = X[ID,K,2,Iz,IF]
    x3 = X[ID,K,3,Iz,IF] 
    xS = SVector{3}(x1, x2 ,x3)
    _,uP,vP,wP = Profile(xS,time)
    lon,lat,_= Grids.cart2sphere(x1,x2,x3)
    u[K,Iz,ind] = uP
    v[K,Iz,ind] = vP
  end
end

@kernel inbounds = true function RhouvFunCKernel!(Profile,Rho,u,v,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[ID,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,uP,vP,_ = Profile(xS,time)
    Rho[Iz,ind] = RhoP
    u[1,Iz,ind] = uP
    v[1,Iz,ind] = vP
  end
end

@kernel inbounds = true function RhoFunCKernel!(Profile,Rho,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ = Profile(xS,time)
    Rho[1,Iz,ind] = RhoP
  end
end

@kernel inbounds = true function ThBGFunCKernel!(Profile,ThBG,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,_,_,_,_,_,_,_,_,ThBGP = Profile(xS,time)
    ThBG[1,Iz,ind] = ThBGP
  end
end

@kernel inbounds = true function RhoFunCDGKernel!(Profile,Rho,time,@Const(Glob),@Const(X),Param,Phys)

  ID, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  M = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz
    ind = Glob[ID,IF]
    x1 = X[ID,K,1,Iz,IF]
    x2 = X[ID,K,2,Iz,IF]
    x3 = X[ID,K,3,Iz,IF]
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ = Profile(xS,time)
    Rho[K,Iz,ind] = RhoP
  end
end

@kernel inbounds = true function TrFunCKernel!(Profile,Tr,time,@Const(Glob),@Const(X),Param,Phys)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,TrP = Profile(xS,time)
    Tr[1,Iz,ind] = RhoP * TrP
  end
end

@kernel inbounds = true function RhoThFunCKernel!(Profile,RhoTh,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,ThP = Profile(xS,time)
    RhoTh[1,Iz,ind] = RhoP * ThP
  end
end

@kernel inbounds = true function RhoThFunCDGKernel!(Profile,RhoTh,time,@Const(Glob),@Const(X))

  I, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = X[I,K,1,Iz,IF]
    x2 = X[I,K,2,Iz,IF]
    x3 = X[I,K,3,Iz,IF]
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,ThP = Profile(xS,time)
    RhoTh[K,Iz,ind] = RhoP * ThP
  end
end

@kernel inbounds = true function ThFunCDGKernel!(Profile,RhoTh,time,@Const(Glob),@Const(X))

  I, K, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = X[I,K,1,Iz,IF]
    x2 = X[I,K,2,Iz,IF]
    x3 = X[I,K,3,Iz,IF]
    xS = SVector{3}(x1, x2 ,x3)
    _,_,_,_ ,ThP = Profile(xS,time)
    RhoTh[K,Iz,ind] = ThP
  end
end

@kernel inbounds = true function RhoEFunCKernel!(Profile,RhoE,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,_,EP = Profile(xS,time)
    RhoE[1,Iz,ind] = RhoP * EP
  end
end

@kernel inbounds = true function RhoIEFunCKernel!(Profile,RhoIE,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_,_,_,IEP = Profile(xS,time)
    RhoIE[1,Iz,ind] = RhoP * IEP
  end
end

@kernel inbounds = true function RhoTFunCKernel!(Profile,RhoV,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_,_,_,_,QvP,QcP = Profile(xS,time)
    RhoV[1,Iz,ind] = RhoP * (QvP + QcP)
  end
end

@kernel inbounds = true function RhoVFunCKernel!(Profile,RhoV,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_,_,_,_,QvP = Profile(xS,time)
    RhoV[1,Iz,ind] = RhoP * QvP
  end
end

@kernel inbounds = true function RhoCFunCKernel!(Profile,RhoC,time,@Const(Glob),@Const(X))

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_,_,_,_,_,QcP = Profile(xS,time)
    RhoC[1,Iz,ind] = RhoP * QcP
  end
end

@kernel inbounds = true function ComputeFunFKernel!(Profile,w,time,@Const(Glob),@Const(X),Param)

  I, iz   = @index(Local, NTuple)
  _,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  if Iz <= Nz - 1 
    ind = Glob[ID,IF]
    x1 = eltype(X)(0.5) * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = eltype(X)(0.5) * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = eltype(X)(0.5) * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,_,_,wP = Profile(xS,time)
    w[1,Iz,ind] = wP
  end
end

#=
@inline function VelSphere2Cart(VelSp,lon,lat)
  VelCa1 = -sin(lon) * VelSp[1] -
    sin(lat) * cos(lon)* VelSp[2] + 
    cos(lat) * cos(lon) * VelSp[3]
  VelCa2 = cos(lon) * VelSp[1] -
    sin(lat) * sin(lon) * VelSp[2] +
    cos(lat) * sin(lon) * VelSp[3]
  VelCa3 =  cos(lat) * VelSp[2] +
    sin(lat) * VelSp[3]
  return SVector{3}(VelCa1, VelCa2 ,VelCa3)
end
=#
