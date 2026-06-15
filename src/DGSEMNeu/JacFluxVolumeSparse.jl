@inline function Index(k,i,j,iz,iF,Pos,M,nz,N,NF)
  ind = k + (iz - 1) * M + (i - 1) * M * nz + (j-1) * M * nz * N + 
    (iF - 1) * M * nz * N * N + (Pos - 1) * M * nz * N * N * NF
end  

@inline function Index(k,ijF,iz,Pos,M,nz,N,NF)
  ind = k + (iz - 1) * M + (ijF - 1) * M * nz + 
    (Pos - 1) * M * nz * N * N * NF
end    

@inline function Insert!(iRow,iCol,v,RowInd,ColInd,Val)
  push!(RowInd,iRow)
  push!(ColInd,iCol)
  push!(Val,v)
end  

function JacFluxVolumeSparseH(DG,Metric,Aux,Phys)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]

  DVT = deepcopy(DG.DVT)
  N = size(DVT,1)
  DVZT = deepcopy(DG.DVZT)
  M = size(DVZT,1)
  DVT = deepcopy(DG.DVT)
  for i = 1 : N
    DVT[i,i] = 0.0  
    DVT[i,i] = sum(DVT[:,i])
  end  
  DVZT = deepcopy(DG.DVZT)
  for i = 1 : M
    DVZT[i,i] = 0.0  
    DVZT[i,i] = sum(DVZT[:,i])
  end  
  dXdxIG = Metric.dXdxI
  nz = size(dXdxIG,5)
  NF = size(dXdxIG,6)
  dXdxI = reshape(dXdxIG,3,3,M,N,N,nz,NF)
  Theta = reshape(Aux[:,:,:,4],M,nz,N,N,NF)
  pRhoTh = reshape(Aux[:,:,:,3],M,nz,N,N,NF)
  Mnz = nz * M
  MnzNN = N * N * M * nz
  ind = 0
  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N  
        ind += 1  
        for iz = 1 : nz  
          for k = 1 : M  
            # Rho RhoTheta 
            iRowRho = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            iRowRhoTh = Index(k,i,j,iz,iF,RhoThPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,uPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,1,k,l,j,iz,iF] * DVT[l,i]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,l,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,l,j,iz,iF,vPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,2,k,l,j,iz,iF] * DVT[l,i]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,l,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,l,j,iz,iF,wPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,3,k,l,j,iz,iF] * DVT[l,i]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,l,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
            end  
#           j direction 
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,uPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,1,k,i,l,iz,iF] * DVT[l,j]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,i,l,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,i,l,iz,iF,vPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,2,k,i,l,iz,iF] * DVT[l,j]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,i,l,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(k,i,l,iz,iF,wPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,3,k,i,l,iz,iF] * DVT[l,j]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[k,iz,i,l,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
            end      
            # u  
            iRow = Index(k,i,j,iz,iF,uPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,1,k,l,j,iz,iF] * pRhoTh[k,iz,l,j,iF] * DVT[l,i]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end  
#           j direction
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,1,k,i,l,iz,iF] * pRhoTh[k,iz,i,l,iF] * DVT[l,j]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end  

            # v  
            iRow = Index(k,i,j,iz,iF,vPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,2,k,l,j,iz,iF] * pRhoTh[k,iz,l,j,iF] * DVT[l,i]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end 
#           j direction
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,2,k,i,l,iz,iF] * pRhoTh[k,iz,i,l,iF] * DVT[l,j]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end 

            # w  
            iRow = Index(k,i,j,iz,iF,wPos,M,nz,N,NF)
#           i direction
            for l = 1 : N
              iCol = Index(k,l,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[1,3,k,l,j,iz,iF] * pRhoTh[k,iz,l,j,iF] * DVT[l,i]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end
#           j direction
            for l = 1 : N
              iCol = Index(k,i,l,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[2,3,k,i,l,iz,iF] * pRhoTh[k,iz,i,l,iF] * DVT[l,j]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end
          end        
        end        
      end        
    end        
  end        
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  return Jac
end

function JacFluxVolumeSparseV(DG,Metric,Aux,Phys)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]

  DVT = deepcopy(DG.DVT)
  N = size(DVT,1)
  DVZT = deepcopy(DG.DVZT)
  M = size(DVZT,1)
  DVT = deepcopy(DG.DVT)
  for i = 1 : N
    DVT[i,i] = 0.0  
    DVT[i,i] = sum(DVT[:,i])
  end  
  DVZT = deepcopy(DG.DVZT)
  for i = 1 : M
    DVZT[i,i] = 0.0  
    DVZT[i,i] = sum(DVZT[:,i])
  end  
  dXdxIG = Metric.dXdxI
  nz = size(dXdxIG,5)
  NF = size(dXdxIG,6)
  dXdxI = reshape(dXdxIG,3,3,M,N,N,nz,NF)
  Theta = reshape(Aux[:,:,:,4],M,nz,N,N,NF)
  pRhoTh = reshape(Aux[:,:,:,3],M,nz,N,N,NF)
  Mnz = nz * M
  MnzNN = N * N * M * nz
  ind = 0
  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N  
        ind += 1  
        for iz = 1 : nz  
          for k = 1 : M  
            # Rho RhoTheta 
            iRowRho = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            iRowRhoTh = Index(k,i,j,iz,iF,RhoThPos,M,nz,N,NF)
#           k direction 
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,uPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,1,l,i,j,iz,iF] * DVZT[l,k]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[l,iz,i,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(l,i,j,iz,iF,vPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,2,l,i,j,iz,iF] * DVZT[l,k]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[l,iz,i,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
              iCol = Index(l,i,j,iz,iF,wPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,3,l,i,j,iz,iF] * DVZT[l,k]
              Insert!(iRowRho,iCol,ValLoc,RowInd,ColInd,Val)
              ValLoc *= Theta[l,iz,i,j,iF]
              Insert!(iRowRhoTh,iCol,ValLoc,RowInd,ColInd,Val)
            end         
            # u  
            iRow = Index(k,i,j,iz,iF,uPos,M,nz,N,NF)
#           k direction
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,1,l,i,j,iz,iF] * pRhoTh[l,iz,i,j,iF] * DVZT[l,k]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end  

            # v  
            iRow = Index(k,i,j,iz,iF,vPos,M,nz,N,NF)
#           k direction
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,2,l,i,j,iz,iF] * pRhoTh[l,iz,i,j,iF] * DVZT[l,k]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end         

            # w  
            iRow = Index(k,i,j,iz,iF,wPos,M,nz,N,NF)
#           k direction
            for l = 1 : M
              iCol = Index(l,i,j,iz,iF,RhoThPos,M,nz,N,NF)
              ValLoc = -0.5 * dXdxI[3,3,l,i,j,iz,iF] * pRhoTh[l,iz,i,j,iF] * DVZT[l,k]
              Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
            end        
          end        
        end        
      end        
    end        
  end        
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  return Jac
end

