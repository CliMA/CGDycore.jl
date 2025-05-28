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

function InitJacDG(DG,nz,J,Param)
  N = (DG.OrdPolyZ + 1) * nz
  dSdS,dSdM = DGVertical.DScalarDMomAc(nz,DG,Param.cS)
  dMdS,dMdM = DGVertical.DMomDScalarAc(nz,DG,Param.cS)

  dSdS = spdiagm(reshape(1.0./J,N)) * dSdS
  dSdM = spdiagm(reshape(1.0./J,N)) * dSdM
  dMdS = spdiagm(reshape(1.0./J,N)) * dMdS
  dMdM = spdiagm(reshape(1.0./J,N)) * dMdM
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
# Jac = [spzeros(N,N) dSdM  spzeros(N,N)
       -Phys.Grav * sparse(I,N,N) dMdM  dMdS * diagm(dpdRhoTh)
       spzeros(N,N) dSdM* diagm(Th)  diagm(Th) * dSdS* diagm(dpdRhoTh)]
end         

function JacAccoustic(U,dSdM,dMdS,fac,Phys)
  N = size(dSdM,1)
  spId = sparse(I,N,N)
  spZero = spzeros(N,N)
  Jac = [spId -fac*dSdM
         -fac*dMdS spId]
end         
