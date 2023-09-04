@kernel function InterpolateKernel!(cCell,@Const(c),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)


  ColumnTilesDim = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  @uniform N = size(Glob,1)

  
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    @inbounds cCell[I,J,Iz,IF] = 0
    for jP = 1 : N
      for iP = 1 : N
        @inbounds ind = Glob[iP,jP,IF]
        @inbounds  cCell[I,J,Iz,IF] += Inter[I,J,iP,jP] * c[Iz,ind]
      end
    end
  end
end

@kernel function InterpolateRhoKernel!(cCell,@Const(c),@Const(RhoC),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)


  ColumnTilesDim = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  @uniform N = size(Glob,1)
  
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    @inbounds cCell[I,J,Iz,IF] = 0
    for jP = 1 : N
      for iP = 1 : N
        @inbounds ind = Glob[iP,jP,IF]
        @inbounds  cCell[I,J,Iz,IF] += Inter[I,J,iP,jP] * c[Iz,ind] / RhoC[Iz,ind]
      end
    end
  end
end

@kernel function InterpolateWBKernel!(cCell,@Const(u),@Const(v),@Const(w),@Const(Inter),@Const(dXdxI),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)


  ColumnTilesDim = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  @uniform N = size(Glob,1)


  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz == 1
    for jP = 1 : N
      for iP = 1 : N
        @inbounds ind = Glob[iP,jP,IF]
        @inbounds w0 = -(u[ind] * dXdxI[3,1,1,iP,jP,Iz,IF] +
          v[ind] * dXdxI[3,2,1,iP,jP,Iz,IF]) / dXdxI[3,3,1,iP,jP,Iz,IF]
        @inbounds  cCell[I,J,Iz,IF] += 0.5 * Inter[I,J,iP,jP] * (w[Iz,ind] + w0)
      end
    end
  elseif Iz <= Nz
    @inbounds cCell[I,J,Iz,IF] = 0
    for jP = 1 : N
      for iP = 1 : N
        @inbounds ind = Glob[iP,jP,IF]
        @inbounds  cCell[I,J,Iz,IF] += 0.5 * Inter[I,J,iP,jP] * (w[Iz,ind] + w[Iz-1,ind])
      end
    end
  end
end

function InterpolateGPU!(cCell,c,Inter,Glob)

  backend = get_backend(c)
  FT = eltype(c)

  OrdPrint = size(Inter,1)
  NF = size(Glob,3)
  Nz = size(c,1)

# Ranges
  NzG = min(div(1024,OrdPrint*OrdPrint),Nz)
  group = (OrdPrint, OrdPrint, NzG, 1)
  ndrange = (OrdPrint, OrdPrint, Nz, NF)

  KInterpolateKernel! = InterpolateKernel!(backend,group)
  KInterpolateKernel!(cCell,c,Inter,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

function InterpolateRhoGPU!(cCell,c,Rho,Inter,Glob)

  backend = get_backend(c)
  FT = eltype(c)

  OrdPrint = size(Inter,1)
  NF = size(Glob,3)
  Nz = size(c,1)

# Ranges
  NzG = min(div(1024,OrdPrint*OrdPrint),Nz)
  group = (OrdPrint, OrdPrint, NzG, 1)
  ndrange = (OrdPrint, OrdPrint, Nz, NF)

  KInterpolateRhoKernel! = InterpolateRhoKernel!(backend,group)
  KInterpolateRhoKernel!(cCell,c,Rho,Inter,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

function InterpolateWBGPU!(cCell,u,v,w,Inter,dXdxI,Glob)

  backend = get_backend(w)
  FT = eltype(w)

  OrdPrint = size(Inter,1)
  NF = size(Glob,3)
  Nz = size(w,1)

# Ranges
  NzG = min(div(1024,OrdPrint*OrdPrint),Nz)
  group = (OrdPrint, OrdPrint, NzG, 1)
  ndrange = (OrdPrint, OrdPrint, Nz, NF)

  KInterpolateWBKernel! = InterpolateWBKernel!(backend,group)
  @views KInterpolateWBKernel!(cCell,u[1,:],v[1,:],w,Inter,dXdxI,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end
