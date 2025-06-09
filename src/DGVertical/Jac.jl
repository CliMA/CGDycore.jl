function DScalarDMomAc(NZ,DG,cS)
  
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
      push!(Val,-0.5/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-0.5/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,0.5/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,0.5/DG.wZ[1])
    end  
  end
  dSdM = sparse(RowInd, ColInd, Val)

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-0.5/cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+0.5/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,0.5/cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-0.5/cS/DG.wZ[1])
    end  
  end
  dSdS = sparse(RowInd, ColInd, Val,N,N)
  return dSdS,dSdM
end

function DMomDScalarAc(NZ,DG,cS)
  
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
      push!(Val,-0.5/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-0.5/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,0.5/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,0.5/DG.wZ[1])
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
  dMdS = sparse(RowInd, ColInd, Val)
  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]
  for iZ = 1 : NZ
    if iZ < NZ
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M)  
      push!(Val,-0.5*cS/DG.wZ[1])
      push!(RowInd,iZ*M)  
      push!(ColInd,iZ*M+1)  
      push!(Val,+0.5*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M)  
      push!(Val,0.5*cS/DG.wZ[1])
      push!(RowInd,iZ*M+1)  
      push!(ColInd,iZ*M+1)  
      push!(Val,-0.5*cS/DG.wZ[1])
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
  @show N,DG.OrdPolyZ + 1,nz
  dSdS,dSdM = DGVertical.DScalarDMomAc(nz,DG,Param.cS)
  dMdS,dMdM = DGVertical.DMomDScalarAc(nz,DG,Param.cS)
  @show size(dSdS),size(dSdM)
  return dSdS,dSdM,dMdS,dMdM
end  

function JacDG(U,dSdS,dSdM,dMdS,dMdM,Phys)
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 3
  Th = reshape(U[:,:,ThPos]./U[:,:,RhoPos],N)
  dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
    (Phys.Rd * U[:,:,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
  Jac = [spzeros(N,N) dSdM  dSdS* diagm(dpdRhoTh)
       -Phys.Grav * sparse(I,N,N) dMdM  dMdS * diagm(dpdRhoTh)
       spzeros(N,N) dSdM* diagm(Th)  diagm(Th) * dSdS* diagm(dpdRhoTh)]
end         

function JacDG1(U,DG,dSdS,dSdM,dMdS,dMdM,J,Phys)
  FTB = eltype(U)
  N = size(dSdM,1)
  RhoPos = 1
  ThPos = 3
  nz = size(U,2)
  NF = size(J,4)
  Jac = Array{SparseMatrixCSC}(undef,1)
  ID = 1  
  @show size(J),N
  @show size(U)
  diagJ = spdiagm(1.0./reshape(J[:,:],N))
  Th = reshape(U[:,:,ThPos]./U[:,:,RhoPos],N)
  dpdRhoTh = reshape( FTB(1) / (FTB(1) - Phys.kappa) * Phys.Rd *
    (Phys.Rd * U[:,:,ThPos] ./ Phys.p0).^(Phys.kappa / (1.0 - Phys.kappa)),N)
  Jac[ID] = [spzeros(N,N)              diagJ * dSdM              diagJ* dSdS * diagm(dpdRhoTh)
            -Phys.Grav * sparse(I,N,N) diagJ * dMdM              diagJ* dMdS * diagm(dpdRhoTh)
#            spzeros(N,N)              diagJ * dSdM * diagm(Th)  diagm(Th) * diagJ * dSdS * diagm(dpdRhoTh)]
             spzeros(N,N)              diagm(Th) * diagJ * dSdM  diagm(Th) * diagJ * dSdS * diagm(dpdRhoTh)]
  return Jac     
end

function JacAccoustic(U,dSdM,dMdS,fac,Phys)
  N = size(dSdM,1)
  spId = sparse(I,N,N)
  spZero = spzeros(N,N)
  Jac = [spId -fac*dSdM
         -fac*dMdS spId]
end         
