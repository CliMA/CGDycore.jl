function JacGravitySparseH(DG,Metric,Aux,Phys)

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
  Geo = reshape(Aux[:,:,:,2],M,nz,N,N,NF)
  Mnz = nz * M
  MnzNN = N * N * M * nz
  ind = 0
  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N  
        ind += 1  
        for iz = 1 : nz  
          for k = 1 : M  
            # u  
            iRow = Index(k,i,j,iz,iF,uPos,M,nz,N,NF)
#           i direction
            ValLoc_i = 0.0  
            for l = 1 : N
              if l != i
                ValLoc = 0.25 * DVT[l,i] * (Geo[k,iz,i,j,iF] - Geo[k,iz,l,j,iF])  
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[1,1,k,l,j,iz,iF]
                iCol = Index(k,l,j,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end    
            end  
            ValLoc_i *= dXdxI[1,1,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)

#           j direction 
            ValLoc_i = 0.0
            for l = 1 : N
              if l != j
                ValLoc = 0.25 * DVT[l,j] * (Geo[k,iz,i,j,iF] - Geo[k,iz,i,l,iF])
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[2,1,k,i,l,iz,iF]
                iCol = Index(k,i,l,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end
            end
            ValLoc_i *= dXdxI[2,1,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)

            # v  
            iRow = Index(k,i,j,iz,iF,vPos,M,nz,N,NF)
#           i direction
            ValLoc_i = 0.0  
            for l = 1 : N
              if l != i
                ValLoc = 0.25 * DVT[l,i] * (Geo[k,iz,i,j,iF] - Geo[k,iz,l,j,iF])  
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[1,2,k,l,j,iz,iF]
                iCol = Index(k,l,j,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end    
            end  
            ValLoc_i *= dXdxI[1,2,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)
#           j direction 
            ValLoc_i = 0.0
            for l = 1 : N
              if l != j
                ValLoc = 0.25 * DVT[l,j] * (Geo[k,iz,i,j,iF] - Geo[k,iz,i,l,iF])
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[2,2,k,i,l,iz,iF]
                iCol = Index(k,i,l,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end
            end
            ValLoc_i *= dXdxI[2,2,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)
            # w  
            iRow = Index(k,i,j,iz,iF,wPos,M,nz,N,NF)
#           i direction
            ValLoc_i = 0.0  
            for l = 1 : N
              if l != i
                ValLoc = 0.25 * DVT[l,i] * (Geo[k,iz,i,j,iF] - Geo[k,iz,l,j,iF])  
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[1,3,k,l,j,iz,iF]
                iCol = Index(k,l,j,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end    
            end  
            ValLoc_i *= dXdxI[1,3,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)
#           j direction 
            ValLoc_i = 0.0
            for l = 1 : N
              if l != j
                ValLoc = 0.25 * DVT[l,j] * (Geo[k,iz,i,j,iF] - Geo[k,iz,i,l,iF])
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[2,3,k,i,l,iz,iF]
                iCol = Index(k,i,l,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end
            end
            ValLoc_i *= dXdxI[2,3,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)
          end        
        end        
      end        
    end        
  end        
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  return Jac
end

function JacGravitySparseV(DG,Metric,Aux,Phys)

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
  Geo = reshape(Aux[:,:,:,2],M,nz,N,N,NF)
  Mnz = nz * M
  MnzNN = N * N * M * nz
  ind = 0
  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N  
        ind += 1  
        for iz = 1 : nz  
          for k = 1 : M  
            # u  
            iRow = Index(k,i,j,iz,iF,uPos,M,nz,N,NF)
#           k direction 
            ValLoc_i = 0.0
            for l = 1 : M
              if l != k
                ValLoc = 0.25 * DVZT[l,k] * (Geo[k,iz,i,j,iF] - Geo[l,iz,i,j,iF])
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[3,1,l,i,j,iz,iF]
                iCol = Index(l,i,j,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end
            end         
            ValLoc_i *= dXdxI[3,1,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)

            # v  
            iRow = Index(k,i,j,iz,iF,vPos,M,nz,N,NF)
#           k direction 
            ValLoc_i = 0.0
            for l = 1 : M
              if l != k
                ValLoc = 0.25 * DVZT[l,k] * (Geo[k,iz,i,j,iF] - Geo[l,iz,i,j,iF])
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[3,2,l,i,j,iz,iF]
                iCol = Index(l,i,j,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end
            end         
            ValLoc_i *= dXdxI[3,2,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)
            # w  
            iRow = Index(k,i,j,iz,iF,wPos,M,nz,N,NF)
#           k direction 
            ValLoc_i = 0.0
            for l = 1 : M
              if l != k
                ValLoc = 0.25 * DVZT[l,k] * (Geo[k,iz,i,j,iF] - Geo[l,iz,i,j,iF])
                ValLoc_i += ValLoc
                ValLoc *= dXdxI[3,3,l,i,j,iz,iF]
                iCol = Index(l,i,j,iz,iF,RhoPos,M,nz,N,NF)
                Insert!(iRow,iCol,ValLoc,RowInd,ColInd,Val)
              end
            end         
            ValLoc_i *= dXdxI[3,3,k,i,j,iz,iF]
            iCol = Index(k,i,j,iz,iF,RhoPos,M,nz,N,NF)
            Insert!(iRow,iCol,ValLoc_i,RowInd,ColInd,Val)
          end        
        end        
      end        
    end        
  end        
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  return Jac
end

