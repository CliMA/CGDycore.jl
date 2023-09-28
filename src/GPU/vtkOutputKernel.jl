@kernel function InterpolateRotKernel!(Rot,@Const(uC),@Const(vC),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(Inter),@Const(Glob))

  gi, gz, gF = @index(Group, NTuple)
  ID, iz   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform size(D,1)
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  RotLoc = @localmem eltype(Rot) (N,N,ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    RotLoc[I,J] = 0
  end
  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz 
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    @inbounds ind = Glob[ID,IF]
    @inbounds tempx = ((dXdxI[1,1,1,ID,Iz,IF] + dXdxI[1,1,2,ID,Iz,IF]) * vC[Iz,ind] - 
      (dXdxI[1,2,1,ID,Iz,IF] + dXdxI[1,2,2,ID,Iz,IF]) * uC[Iz,ind])
    @inbounds tempy = ((dXdxI[2,1,1,ID,Iz,IF] + dXdxI[2,1,2,ID,Iz,IF]) * vC[Iz,ind] - 
      (dXdxI[2,2,1,ID,Iz,IF] + dXdxI[2,2,2,ID,Iz,IF]) * uC[Iz,ind])
    for k = 1 : N
      @inbounds @atomic RotLoc[k,J,iz] += D[k,I] * tempx
      @inbounds @atomic RotLoc[I,k,iz] += D[k,J] * tempy
    end
  end
  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz 
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    @inbounds temp = RotLoc[I,J,iz] / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    for j = 1 : size(Inter,1) 
      for i = 1 : size(Inter,2) 
        @inbounds @atomic Rot[i,j,Iz,IF] += Inter[i,j,I,J] * temp
      end    
    end    
  end    

end  

@kernel function InterpolateKernel!(cCell,@Const(c),@Const(Inter),@Const(Glob),
  ::Val{BANK}=Val(1)) where BANK
  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)


  ColumnTilesDim = @uniform @groupsize()[3]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  @uniform N = size(Inter,3)

  
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    @inbounds cCell[I,J,Iz,IF] = 0
    iD = 0
    for jP = 1 : N
      for iP = 1 : N
        iD += 1  
        @inbounds ind = Glob[iD,IF]
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
  @uniform N = size(Inter,3)
  
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    @inbounds cCell[I,J,Iz,IF] = 0
    ID = 0
    for jP = 1 : N
      for iP = 1 : N
        ID += 1
        @inbounds ind = Glob[ID,IF]
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

function InterpolateVortGPU!(Vort,U,Inter,FE,Metric)

  backend = get_backend(U)
  FT = eltype(Vort)

  OrdPrint = size(Inter,1)
  DoF = FE.DoF
  Glob = FE.Glob
  NF = size(Glob,2)
  Nz = size(U,1)

# Ranges
  NzG = min(div(1024,DoF),Nz)
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


function InterpolateGPU!(cCell,c,Inter,Glob)

  backend = get_backend(c)
  FT = eltype(c)

  OrdPrint = size(Inter,1)
  NF = size(Glob,2)
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
  NF = size(Glob,2)
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
