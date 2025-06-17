function DScalarDMomAc(NZ,DG,cS)
  
  fac = 0.5
  M = DG.OrdPolyZ + 1
  N = NZ * M
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    for i = 1 : M
      for j = 1 : M
        push!(RowInd,i+(iZ-1)*M)
        push!(ColInd,j+(iZ-1)*M)
        push!(Val,-DG.DWZ[i,j])
      end
    end  
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,fac/DG.wZ[1])
    end  
  end
  dSdM = sparse(RowInd, ColInd, Val,N,N)

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/cS/DG.wZ[1])
    end  
  end
  dSdS = sparse(RowInd, ColInd, Val,N,N)
  return dSdS,dSdM
end

function DMomDScalarAc(NZ,DG,cS)
  
  fac = 0.5
  M = DG.OrdPolyZ + 1
  N = NZ * M
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    for i = 1 : M
      for j = 1 : M
        push!(RowInd,i+(iZ-1)*M)
        push!(ColInd,j+(iZ-1)*M)
        push!(Val,-DG.DWZ[i,j])
      end
    end  
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,fac/DG.wZ[1])
    end  
    if iZ == 1
      push!(RowInd,1)  
      push!(ColInd,1)  
      push!(Val,1/DG.wZ[1])
    end    
    if iZ == NZ
      push!(RowInd,M*NZ)  
      push!(ColInd,M*NZ)  
      push!(Val,-1/DG.wZ[1])
    end    
  end
  dMdS = sparse(RowInd, ColInd, Val,N,N)
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,fac*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-fac*cS/DG.wZ[1])
    end  
    if iZ == 1
      push!(RowInd,1)  
      push!(ColInd,1)  
      push!(Val,-cS/DG.wZ[1])
    end    
    if iZ == NZ
      push!(RowInd,M*NZ)  
      push!(ColInd,M*NZ)  
      push!(Val,-cS/DG.wZ[1])
    end    
  end  
  dMdM = sparse(RowInd, ColInd, Val,N,N)
  return dMdS,dMdM
end

function InitJacDG(DG,nz,Param)
  N = (DG.OrdPolyZ + 1) * nz
  dSdS,dSdM = DScalarDMomAc(nz,DG,Param.cS)
  dMdS,dMdM = DMomDScalarAc(nz,DG,Param.cS)
  return dSdS,dSdM,dMdS,dMdM
end  

function JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,z,Phys)
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 5
  nz = size(U,2)
  M = size(U,1)
  oneM = ones(M)
  NF = size(z,3)
  JacLU = Array{SparseArrays.UMFPACK.UmfpackLU}(undef,size(U,3))
  for ID = 1 : DG.NumI
    @views zCol = z[:,ID]
    diagz = spdiagm(2.0 ./ reshape(vec(oneM*zCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Jac = [sparse(fac*I,N,N) -diagz * dSdM              -diagz* dSdS * diagm(dpdRhoTh)
           sparse((0.5 * Phys.Grav)*I,N,N) sparse(fac*I,N,N) - diagz * dMdM -diagz* dMdS * diagm(dpdRhoTh)
           spzeros(N,N) -diagz * dSdM * diagm(Th)  sparse(fac*I,N,N) - diagz * diagm(Th) * dSdS * diagm(dpdRhoTh)]
    JacLU[ID] = lu(Jac)           
  end
  return JacLU
end

function JacDGT(U,DG,fac,dSdS,dSdM,dMdS,dMdM,z,Phys)
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 5
  nz = size(U,2)
  M = size(U,1)
  oneM = ones(M)
  NF = size(z,3)
  JacLU = Array{SparseArrays.UMFPACK.UmfpackLU}(undef,size(U,3))
    ID = 1  
    @views zCol = z[:,ID]
    diagz = spdiagm(2.0 ./ reshape(vec(oneM*zCol'),N))
    Th = reshape(U[:,:,ID,ThPos]./U[:,:,ID,RhoPos],N)
    dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
      (Phys.Rd * U[:,:,ID,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
    Jac = [sparse(fac*I,N,N)  diagz* dSdS * diagm(dpdRhoTh) -diagz * dSdM
           spzeros(N,N) sparse(fac*I,N,N)-diagz*diagm(Th)*dSdS*diagm(dpdRhoTh)  -diagz*dSdM*diagm(Th)
           sparse((0.5 * Phys.Grav)*I,N,N) -diagz*dMdS*diagm(dpdRhoTh) sparse(fac*I,N,N)-diagz*dMdM]
    JacLU[ID] = lu(Jac)
  return JacLU,Jac
end
