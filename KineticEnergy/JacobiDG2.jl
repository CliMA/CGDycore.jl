function JacobiDG2(xP,zP,Fe)
  ksi=Fe.xw
  zeta=Fe.zw
  nX=size(ksi,1)
  nZ=size(zeta,1)
  X=zeros(nX,nZ,2)
  dXdx=zeros(nX,nZ,2,2)
  dXdxI=zeros(nX,nZ,2,2)
  norm_nSB=zeros(nX)
  nSB=zeros(nX,2)
  VolSurfB=zeros(nX)
  norm_nST=zeros(nX)
  nST=zeros(nX,2)
  VolSurfT=zeros(nX)
  J=zeros(nX,nZ)
  for j=1:nZ
    for i=1:nX
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

  for i=1:nX
    for j=1:nZ
      J[i,j]=det(reshape(dXdx[i,j,:,:],2,2))
      dXdxI[i,j,:,:]=inv(reshape(dXdx[i,j,:,:],2,2))*J[i,j]
    end
  end
  for i=1:nX
    norm_nSB[i] = sqrt(dXdxI[i,1,2,1] * dXdxI[i,1,2,1] + dXdxI[i,1,2,2] * dXdxI[i,1,2,2]) 
    nSB[i,1] = dXdxI[i,1,2,1] / norm_nSB[i]
    nSB[i,2] = dXdxI[i,1,2,2] / norm_nSB[i]
    VolSurfB[i] = sqrt(dXdx[i,1,1,1] * dXdx[i,1,1,1] + dXdx[i,1,2,1] * dXdx[i,1,2,1]) 
    norm_nST[i] = sqrt(dXdxI[i,nZ,2,1] * dXdxI[i,nZ,2,1] + dXdxI[i,nZ,2,2] * dXdxI[i,nZ,2,2]) 
    nST[i,1] = dXdxI[i,nZ,2,1] / norm_nST[i]
    nST[i,2] = dXdxI[i,nZ,2,2] / norm_nST[i]
    VolSurfT[i] = sqrt(dXdx[i,nZ,1,1] * dXdx[i,nZ,1,1] + dXdx[i,1,2,1] * dXdx[i,nZ,2,1]) 
  end  
  return (X,J,dXdx,dXdxI,norm_nSB,nSB,VolSurfB,norm_nST,nST,VolSurfT)
end

