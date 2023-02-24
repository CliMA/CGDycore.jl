function JacobiDG2(xP,zP,Fe)
  ksi=Fe.xw
  zeta=Fe.zw
  nx=size(ksi,1)
  nz=size(zeta,1)
  X=zeros(nx,nz,2)
  dXdx=zeros(nx,nz,2,2)
  dXdxI=zeros(nx,nz,2,2)
  norm_nSB=zeros(nx)
  nSB=zeros(nx,2)
  VolSurfB=zeros(nx)
  norm_nST=zeros(nx)
  nST=zeros(nx,2)
  VolSurfT=zeros(nx)
  J=zeros(nx,nz)
  for j=1:nz
    for i=1:nx
      X[i,j,1]=xP[i,j]
      X[i,j,2]=zP[i,j]
    end
  end
  X[:,:,1] = Fe.IntXE2F * xP * Fe.IntZE2FT
  X[:,:,2] = Fe.IntXE2F * zP * Fe.IntZE2FT

  dXdx[:,:,1,1]=Fe.DX*X[:,:,1]
  dXdx[:,:,2,1]=Fe.DX*X[:,:,2]
  dXdx[:,:,1,2]=X[:,:,1]*Fe.DZT
  dXdx[:,:,2,2]=X[:,:,2]*Fe.DZT

  for i=1:nx
    for j=1:nz
      J[i,j]=det(reshape(dXdx[i,j,:,:],2,2))
      dXdxI[i,j,:,:]=inv(reshape(dXdx[i,j,:,:],2,2))*J[i,j]
    end
  end
  for i=1:nx
    norm_nSB[i] = sqrt(dXdxI[i,1,2,1] * dXdxI[i,1,2,1] + dXdxI[i,1,2,2] * dXdxI[i,1,2,2]) 
    nSB[i,1] = dXdxI[i,1,2,1] / norm_nSB[i]
    nSB[i,2] = dXdxI[i,1,2,2] / norm_nSB[i]
    VolSurfB[i] = sqrt(dXdx[i,1,1,1] * dXdx[i,1,1,1] + dXdx[i,1,2,1] * dXdx[i,1,2,1]) 
    norm_nST[i] = sqrt(dXdxI[i,nz,2,1] * dXdxI[i,nz,2,1] + dXdxI[i,nz,2,2] * dXdxI[i,nz,2,2]) 
    nST[i,1] = dXdxI[i,nz,2,1] / norm_nST[i]
    nST[i,2] = dXdxI[i,nz,2,2] / norm_nST[i]
    VolSurfT[i] = sqrt(dXdx[i,nz,1,1] * dXdx[i,nz,1,1] + dXdx[i,1,2,1] * dXdx[i,nz,2,1]) 
  end  
  return (X,J,dXdx,dXdxI,norm_nSB,nSB,VolSurfB,norm_nST,nST,VolSurfT)
end

function JacobiDG3(xP,yP,zP,Fe)
  ksi=Fe.xw
  eta=Fe.yw
  zeta=Fe.zw
  nx=size(ksi,1)
  ny=size(eta,1)
  nz=size(zeta,1)
  X=zeros(nx,ny,nz,3)
  XE=zeros(nx,ny,nz,3)
  dXdx=zeros(nx,ny,nz,3,3)
  dXdxI=zeros(nx,ny,nz,3,3)
  norm_nSB=zeros(nx,ny)
  nSB=zeros(nx,ny,3)
  VolSurfB=zeros(nx,ny)
  norm_nST=zeros(nx,ny)
  nST=zeros(nx,ny,3)
  VolSurfT=zeros(nx,ny)
  J=zeros(nx,ny,nz)
  for i=1:nx
    for j=1:ny  
      for k=1:nz
        XE[i,j,k,1]=xP[i,j,k]
        XE[i,j,k,2]=yP[i,j,k]
        XE[i,j,k,3]=zP[i,j,k]
      end
    end
  end
  ChangeBasis!(X,XE,Fe)

  for k = 1 : nz  
    for j = 1 : ny
      for i = 1 : nx  
        for l = 1 : nx
          @views @. dXdx[i,j,k,:,1] = dXdx[i,j,k,:,1] + Fe.DX[i,l] * X[l,j,k,:]  
        end  
      end
    end
  end  

  for k = 1 : nz  
    for i = 1 : nx  
      for j = 1 : ny
        for l = 1 : ny
          @views @. dXdx[i,j,k,:,2] = dXdx[i,j,k,:,2] + Fe.DY[j,l] * X[i,l,k,:]  
        end  
      end
    end
  end  
  for j = 1 : ny
    for i = 1 : nx  
      for k = 1 : nz  
        for l = 1 : nz
          @views @. dXdx[i,j,k,:,3] = dXdx[i,j,k,:,3] + Fe.DZ[k,l] * X[i,j,l,:]
        end
      end
    end
  end
          
  for i = 1 : nx
    for j = 1 : ny  
      for k = 1 : nz
        J[i,j,k] = det(reshape(dXdx[i,j,k,:,:],3,3))
        dXdxI[i,j,k,:,:] = inv(reshape(dXdx[i,j,k,:,:],3,3))*J[i,j,k]
      end
    end
  end
  for i = 1 : nx
    for j = 1 : ny
#     norm_nSB[i,j] = sqrt(dXdxI[i,1,2,1] * dXdxI[i,1,2,1] + dXdxI[i,1,2,2] * dXdxI[i,1,2,2]) 
#     nSB[i,1] = dXdxI[i,1,2,1] / norm_nSB[i]
#     nSB[i,2] = dXdxI[i,1,2,2] / norm_nSB[i]
#     VolSurfB[i] = sqrt(dXdx[i,1,1,1] * dXdx[i,1,1,1] + dXdx[i,1,2,1] * dXdx[i,1,2,1]) 
#     norm_nST[i] = sqrt(dXdxI[i,nz,2,1] * dXdxI[i,nz,2,1] + dXdxI[i,nz,2,2] * dXdxI[i,nz,2,2]) 
#     nST[i,1] = dXdxI[i,nz,2,1] / norm_nST[i]
#     nST[i,2] = dXdxI[i,nz,2,2] / norm_nST[i]
#     VolSurfT[i] = sqrt(dXdx[i,nz,1,1] * dXdx[i,nz,1,1] + dXdx[i,1,2,1] * dXdx[i,nz,2,1]) 
    end  
  end  
  return (X,J,dXdx,dXdxI,norm_nSB,nSB,VolSurfB,norm_nST,nST,VolSurfT)
end

