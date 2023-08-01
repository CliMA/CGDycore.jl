using CGDycore
using MPI
using Base
using Metal

function DivRhoTrColumnE1!(FRhoTrC,uC,vC,w,RhoTrC,D,dXdxI)

  N = size(uC,1)
  Nz = size(uC,3)

  @. FRhoTrC = 0
  @inbounds for iz = 1 : Nz
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N
        temp1 = -RhoTrC[i,j,iz] * (dXdxI[i,j,1,iz,3,1] * uC[i,j,iz] +
          dXdxI[i,j,1,iz,3,2] * vC[i,j,iz] + dXdxI[i,j,1,iz,3,3] * w[i,j,iz])
        temp2 = -RhoTrC[i,j,iz] * (dXdxI[i,j,2,iz,3,1] * uC[i,j,iz] +
          dXdxI[i,j,2,iz,3,2] * vC[i,j,iz] + dXdxI[i,j,2,iz,3,3] * w[i,j,iz+1])
        FRhoTrC[i,j,iz] += (temp2 - temp1)
        if iz > 1
          FRhoTrC[i,j,iz-1] -= 1/2 * temp1
          FRhoTrC[i,j,iz] += 1/2 * temp2
        end
        if iz < Nz
          FRhoTrC[i,j,iz] -= 1/2 * temp2
          FRhoTrC[i,j,iz+1] += 1/2 * temp1
        end
      end
    end
  end
  @inbounds for iz = 1 : Nz
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N
        tempx = -RhoTrC[i,j,iz] * ((dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * uC[i,j,iz] +
          (dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * vC[i,j,iz] +
          dXdxI[i,j,1,iz,1,3] * w[i,j,iz] + dXdxI[i,j,2,iz,1,3] * w[i,j,iz+1])
        tempy = -RhoTrC[i,j,iz] * ((dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * uC[i,j,iz] +
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * vC[i,j,iz] +
          dXdxI[i,j,1,iz,2,3] * w[i,j,iz] + dXdxI[i,j,2,iz,2,3] * w[i,j,iz+1])
        for k = 1 : N
          FRhoTrC[k,j,iz] += D[k,i] * tempx
          FRhoTrC[i,k,iz] += D[k,j] * tempy
        end
      end
    end
  end
end

function DivRhoTrColumnE2!(FRhoTrC,uC,vC,w,RhoTrC,D,dXdxI)

  Nz = size(F,3)
  N = size(F,1)
  NF = size(F,4)

  zThreads = 30
  threads = (N, N, zThreads)
  zBlocks = cld(N*N*Nz*NF,N*N*zThreads)
  grid  = (1, 1, zBlocks)
  fill(FRhoTrC,0)
  @show grid,threads
# @cuda grid=grid threads=threads DivRhoTrColumnKernel!(FRhoTrC,uC,vC,w,RhoTrC,
#   D,dXdxI,N,Nz,NF)
# Metal.@time @cuda grid=grid threads=threads DivRhoTrColumnKernel!(FRhoTrC,uC,vC,
#   w,RhoTrC,D,dXdxI,N,Nz,NF)

end

function DivRhoTrColumnKernel!(FRhoTrC,uC,vC,w,RhoTrC,D,dXdxI,N,Nz,NF)

#index = thread_position_in_grid_1d()
#  stride = threads_per_threadgroup_1d()

  i = (thread_position_in_grid_1d() - 1) * threads_per_threadgroup_1d() + thread_position_in_threadgroup_1d
  j = (thread_position_in_grid_2d() - 1) * threads_per_threadgroup_1d() + thread_position_in_threadgroup_2d
  l = (thread_position_in_grid_3d() - 1) * threads_per_threadgroup_3d() + thread_position_in_threadgroup_3d

  iz = mod(l,Nz)
  if iz == 0
    iz = Nz
  end
  iF = div(l-iz,Nz) + 1
  
  @inbounds if l <= Nz * NF
    temp1 = -RhoTrC[i,j,iz,iF] * (dXdxI[i,j,1,iz,iF,3,1] * uC[i,j,iz,iF] +
      dXdxI[i,j,1,iz,iF,3,2] * vC[i,j,iz,iF] + dXdxI[i,j,1,iz,iF,3,3] * w[i,j,iz,iF])
    temp2 = -RhoTrC[i,j,iz,iF] * (dXdxI[i,j,2,iz,iF,3,1] * uC[i,j,iz,iF] +
      dXdxI[i,j,2,iz,iF,3,2] * vC[i,j,iz,iF] + dXdxI[i,j,2,iz,iF,3,3] * w[i,j,iz+1,iF])
#   Metal.@atomic FRhoTrC[i,j,iz,iF] += (temp2 - temp1)
    FRhoTrC[i,j,iz,iF] += (temp2 - temp1)
    if iz > 1
#     Metal.@atomic FRhoTrC[i,j,iz-1,iF] -= 1/2 * temp1
#     Metal.@atomic FRhoTrC[i,j,iz,iF] += 1/2 * temp2
      FRhoTrC[i,j,iz-1,iF] -= 1/2 * temp1
      FRhoTrC[i,j,iz,iF] += 1/2 * temp2
    end
    if iz < Nz
#     Metal.@atomic FRhoTrC[i,j,iz,iF] -= 1/2 * temp2
#     Metal.@atomic FRhoTrC[i,j,iz+1,iF] += 1/2 * temp1
      FRhoTrC[i,j,iz,iF] -= 1/2 * temp2
      FRhoTrC[i,j,iz+1,iF] += 1/2 * temp1
    end
    tempx = -RhoTrC[i,j,iz,iF] * ((dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * uC[i,j,iz,iF] +
      (dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * vC[i,j,iz,iF] +
      dXdxI[i,j,1,iz,iF,1,3] * w[i,j,iz,iF] + dXdxI[i,j,2,iz,iF,1,3] * w[i,j,iz+1,iF])
    tempy = -RhoTrC[i,j,iz,iF] * ((dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * uC[i,j,iz,iF] +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * vC[i,j,iz,iF] +
      dXdxI[i,j,1,iz,iF,2,3] * w[i,j,iz,iF] + dXdxI[i,j,2,iz,iF,2,3] * w[i,j,iz+1,iF])
    for k = 1 : N
#     Metal.@atomic FRhoTrC[k,j,iz,iF] += D[k,i] * tempx
#     Metal.@atomic FRhoTrC[i,k,iz,iF] += D[k,j] * tempy
      FRhoTrC[k,j,iz,iF] += D[k,i] * tempx
      FRhoTrC[i,k,iz,iF] += D[k,j] * tempy
    end
  end
  return nothing
end

function DivRhoGradE1!(F,cC,RhoC,D,DW,dXdxI,J)

  Nz = size(F,3)
  N = size(F,1)

  @inbounds for iz = 1 : Nz
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N
        Dxc = 0
        Dyc = 0
        @inbounds for k = 1 : N
          Dxc += D[i,k] * (cC[k,j,iz] / RhoC[k,j,iz])
          Dyc += D[j,k] * (cC[i,k,iz] / RhoC[i,k,iz])
        end
        GradDx = ((dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * Dxc +
          (dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * Dyc) / (J[i,j,1,iz] + J[i,j,2,iz])
        GradDy = ((dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * Dxc +
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * Dyc) / (J[i,j,1,iz] + J[i,j,2,iz])
        tempx = (dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * GradDx +
          (dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * GradDy
        tempy = (dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * GradDx +
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * GradDy
        @inbounds for k = 1 : N
          F[k,j,iz] += DW[k,i] * tempx
          F[i,k,iz] += DW[k,j] * tempy
        end
      end
    end
  end
end

function memset_kernel(array, value)
  i = thread_position_in_grid_1d()
  if i <= length(array)
    @inbounds array[i] = value
  end
  return
end


function DivRhoGradE2!(F,cC,RhoC,DD,DW,dXdxI,J)

  Nz = size(F,3)
  N = size(F,1)
  NF = size(F,4)
  zThreads = 40
  threads = (N, N, zThreads)
  zBlocks = cld(N*N*Nz*NF,N*N*zThreads) 
  grid  = (1, 1, zBlocks)
  fill(F,0)
  @show grid,threads
  a = MtlArray{Float32}(undef, 512)
  @show "Start"
  Metal.@metal threads=512 memset_kernel(a, 42)
  @show "End"
# @metal threads=512 grid=2 DivRhoGradKernel!(F,cC,RhoC,DD,DW,dXdxI,J,N,Nz,NF)
# Metal.@time @cuda grid=grid threads=threads DivRhoGradKernel!(F,cC,RhoC,DD,DW,dXdxI,J,N,Nz,NF)
end

function DivRhoGradKernel!(F,cC,RhoC,D,DW,dXdxI,J,N,Nz,NF)

  i = thread_position_in_threadgroup_1d()
# i = (threadgroup_position_in_grid_3d().x - 1) * threads_per_threadgroup_3d().x + thread_position_in_threadgroup_3d().x
# j = (threadgroup_position_in_grid_3d().y - 1) * threads_per_threadgroup_3d().y + thread_position_in_threadgroup_3d().y
# l = (threadgroup_position_in_grid_3d().z - 1) * threads_per_threadgroup_3d().z + thread_position_in_threadgroup_3d().z
#=
  iz = mod(l,Nz)
  if iz == 0
    iz = Nz
  end
  iF = div(l-iz,Nz) + 1
  
  @inbounds if l <= Nz * NF
    Dxc = 0
    Dyc = 0
    for k = 1 : N
      Dxc = Dxc + D[i,k] * (cC[k,j,iz,iF] / RhoC[k,j,iz,iF])
      Dyc = Dyc + D[j,k] * (cC[i,k,iz,iF] / RhoC[i,k,iz,iF])
    end
    GradDx = ((dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * Dxc +
      (dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * Dyc) / (J[i,j,1,iz,iF] + J[i,j,2,iz,iF])
    GradDy = ((dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * Dxc +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * Dyc) / (J[i,j,1,iz,iF] + J[i,j,2,iz,iF])
    tempx = (dXdxI[i,j,1,iz,iF,1,1] + dXdxI[i,j,2,iz,iF,1,1]) * GradDx +
      (dXdxI[i,j,1,iz,iF,1,2] + dXdxI[i,j,2,iz,iF,1,2]) * GradDy
    tempy = (dXdxI[i,j,1,iz,iF,2,1] + dXdxI[i,j,2,iz,iF,2,1]) * GradDx +
      (dXdxI[i,j,1,iz,iF,2,2] + dXdxI[i,j,2,iz,iF,2,2]) * GradDy
    for k = 1 : N
#     Metal.@atomic F[k,j,iz,iF] = F[k,j,iz,iF] + DW[k,i] * tempx
#     Metal.@atomic F[i,k,iz,iF] = F[i,k,iz,iF] + DW[k,j] * tempy
      F[k,j,iz,iF] = F[k,j,iz,iF] + DW[k,i] * tempx
      F[i,k,iz,iF] = F[i,k,iz,iF] + DW[k,j] * tempy
    end
  end
=#  
  return nothing
end


Nz = 60
NF = 100
OrdPoly = 4

(DW,DS,) = CGDycore.DerivativeMatrixSingle(OrdPoly)
DW32 = zeros(Float32,OrdPoly+1,OrdPoly+1)
DS32 = zeros(Float32,OrdPoly+1,OrdPoly+1)
@. DW32 = Float32(DW)
@. DS32 = Float32(DS)
F = zeros(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
cC = rand(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
uC = rand(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
vC = rand(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
w = rand(Float32,OrdPoly+1,OrdPoly+1,Nz+1,NF)
RhoC = ones(Float32,OrdPoly+1,OrdPoly+1,Nz,NF)
dXdxI = rand(Float32,OrdPoly+1,OrdPoly+1,2,Nz,NF,3,3)
J = ones(Float32,OrdPoly+1,OrdPoly+1,2,Nz,NF)

dDW = MtlArray(DW32)
copyto!(dDW,DW32)
dDS = MtlArray(DS32)
copyto!(dDS,DS32)
dF = MtlArray(F)
copyto!(dF,F)
dcC = MtlArray(cC)
copyto!(dcC,cC)
duC = MtlArray(uC)
copyto!(duC,uC)
dvC = MtlArray(vC)
copyto!(dvC,vC)
dw = MtlArray(w)
copyto!(dw,w)
dRhoC = MtlArray(RhoC)
copyto!(dRhoC,RhoC)
ddXdxI = MtlArray(dXdxI)
copyto!(ddXdxI,dXdxI)
dJ = MtlArray(J)
copyto!(dJ,J)


@. F = 0
for iF = 1 : NF
  @views DivRhoGradE1!(F[:,:,:,iF],cC[:,:,:,iF],RhoC[:,:,:,iF],DS32,DW32,
    dXdxI[:,:,:,:,iF,:,:],J[:,:,:,:,iF])
end
@time for iF = 1 : NF
  @views DivRhoGradE1!(F[:,:,:,iF],cC[:,:,:,iF],RhoC[:,:,:,iF],DS32,DW32,
    dXdxI[:,:,:,:,iF,:,:],J[:,:,:,:,iF])
end
@show "CPU",sum(abs.(F))
DivRhoGradE2!(dF,dcC,dRhoC,dDS,dDW,ddXdxI,dJ)
copyto!(F,dF)
@show "GPU",sum(abs.(F));

@. F = 0
for iF = 1 : NF
  @views DivRhoTrColumnE1!(F[:,:,:,iF],uC[:,:,:,iF],vC[:,:,:,iF],w[:,:,:,iF],
    cC[:,:,:,iF],DS32,dXdxI[:,:,:,:,iF,:,:])
end
@time for iF = 1 : NF
  @views DivRhoTrColumnE1!(F[:,:,:,iF],uC[:,:,:,iF],vC[:,:,:,iF],w[:,:,:,iF],
    cC[:,:,:,iF],DS32,dXdxI[:,:,:,:,iF,:,:])
end
@show "CPU",sum(abs.(F))
DivRhoTrColumnE2!(dF,duC,dvC,dw,dcC,dDS,ddXdxI)
copyto!(F,dF)
@show "GPU",sum(abs.(F));
