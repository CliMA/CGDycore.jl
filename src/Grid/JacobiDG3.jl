function JacobiDG3(DG,F,z,Topo,Topography,zs)
ksi=DG.xw;
eta=DG.xw;
zeta=DG.xwZ;
z1=z[1];
z2=z[2];
n=DG.OrdPoly+1;
nz=DG.OrdPolyZ+1;
X=zeros(n,n,nz,3);
dXdx=zeros(n,n,nz,3,3);
dXdxI=zeros(n,n,nz,3,3);
J=zeros(n,n,nz);
theta=zeros(n,n);
for j=1:n
  for i=1:n
    for k=1:nz
    X[i,j,k,1]=0.25*((1-ksi[i])*(1-eta[j])*F.P[1].x+
                    (1+ksi[i])*(1-eta[j])*F.P[2].x+
                    (1+ksi[i])*(1+eta[j])*F.P[3].x+
                    (1-ksi[i])*(1+eta[j])*F.P[4].x);
    X[i,j,k,2]=0.25*((1-ksi[i])*(1-eta[j])*F.P[1].y+
                    (1+ksi[i])*(1-eta[j])*F.P[2].y+
                    (1+ksi[i])*(1+eta[j])*F.P[3].y+
                    (1-ksi[i])*(1+eta[j])*F.P[4].y);
    z=0.5*((1-zeta[k])*z1+(1+zeta[k])*z2);
    (X[i,j,k,3],dXdx[i,j,k,3,3])=Topo(X[i,j,k,1],X[i,j,k,2],z,z,Topography,zs);
    dXdx[i,j,k,3,3]=0.5*dXdx[i,j,k,3,3]*(z2-z1);
    end
  end
end

for k=1:nz
  dXdx[:,:,k,1,1]=DG.DS*X[:,:,k,1];
  dXdx[:,:,k,2,1]=DG.DS*X[:,:,k,2];
  dXdx[:,:,k,3,1]=DG.DS*X[:,:,k,3];
  dXdx[:,:,k,1,2]=reshape(X[:,:,k,1],n,n)*DG.DST;
  dXdx[:,:,k,2,2]=reshape(X[:,:,k,2],n,n)*DG.DST;
  dXdx[:,:,k,3,2]=reshape(X[:,:,k,3],n,n)*DG.DST;
end
for j=1:n
  for i=1:n
    for k=1:nz
      J[i,j,k]=det(reshape(dXdx[i,j,k,:,:],3,3));
      dXdxI[i,j,k,:,:]=inv(reshape(dXdx[i,j,k,:,:],3,3))*J[i,j,k];
    end
  end
end
return X,J,dXdx,dXdxI,theta
end

function JacobiDG3Neu(DG,F,z,Topo,Topography,zs)
ksi=DG.xe;
eta=DG.xe;
zeta=DG.xwZ;
z1=z[1];
z2=z[2];
n=DG.OrdPoly+1;
nz=DG.OrdPolyZ+1;
XE=zeros(n,n,nz,3);
X=zeros(n,n,nz,3);
dXdx=zeros(n,n,nz,3,3);
dXdxI=zeros(n,n,nz,3,3);
J=zeros(n,n,nz);
theta=zeros(n,n);
for j=1:n
  for i=1:n
    for k=1:nz
    XE[i,j,k,1]=0.25*((1-ksi[i])*(1-eta[j])*F.P[1].x+
                    (1+ksi[i])*(1-eta[j])*F.P[2].x+
                    (1+ksi[i])*(1+eta[j])*F.P[3].x+
                    (1-ksi[i])*(1+eta[j])*F.P[4].x);
    XE[i,j,k,2]=0.25*((1-ksi[i])*(1-eta[j])*F.P[1].y+
                    (1+ksi[i])*(1-eta[j])*F.P[2].y+
                    (1+ksi[i])*(1+eta[j])*F.P[3].y+
                    (1-ksi[i])*(1+eta[j])*F.P[4].y);
    z=0.5*((1-zeta[k])*z1+(1+zeta[k])*z2);
    (XE[i,j,k,3],dXdx[i,j,k,3,3])=Topo(XE[i,j,k,1],XE[i,j,k,2],z,z,Topography,zs);
    dXdx[i,j,k,3,3]=0.5*dXdx[i,j,k,3,3]*(z2-z1);
    dXdx[i,j,k,3,3]=0.0
    end
  end
end

  ChangeBasis!(X,XE,DG)
  for k = 1 : nz
    for j = 1 : n
      for i = 1 : n
        for l = 1 : n
          @views @. dXdx[i,j,k,:,1] = dXdx[i,j,k,:,1] + DG.DS[i,l] * X[l,j,k,:]
        end
      end
    end
  end

  for k = 1 : nz
    for i = 1 : n
      for j = 1 : n
        for l = 1 : n
          @views @. dXdx[i,j,k,:,2] = dXdx[i,j,k,:,2] + DG.DS[j,l] * X[i,l,k,:]
        end
      end
    end
  end
  for j = 1 : n
    for i = 1 : n
      for k = 1 : nz
        for l = 1 : nz
          @views @. dXdx[i,j,k,:,3] = dXdx[i,j,k,:,3] + DG.DSZ[k,l] * X[i,j,l,:]
        end
      end
    end
  end

for j=1:n
  for i=1:n
    for k=1:nz
      J[i,j,k]=det(reshape(dXdx[i,j,k,:,:],3,3));
      dXdxI[i,j,k,:,:]=inv(reshape(dXdx[i,j,k,:,:],3,3))*J[i,j,k];
    end
  end
end
return X,J,dXdx,dXdxI,theta
end


function ChangeBasis!(XOut,XIn,DG)

  nxOut = size(XOut,1)
  nyOut = size(XOut,2)
  nzOut = size(XOut,3)
  nxIn = size(XIn,1)
  nyIn = size(XIn,2)
  nzIn = size(XIn,3)

  Buf1 = zeros(nxOut,nyIn,nzIn,3)
  Buf2 = zeros(nxOut,nyOut,nzIn,3)

  for kIn = 1 : nzIn
    for jIn = 1 : nyIn
      for iIn = 1 : nxIn
        for iOut = 1 : nxOut
          @views @.  Buf1[iOut,jIn,kIn,:] = Buf1[iOut,jIn,kIn,:] + 
            DG.IntXE2F[iOut,iIn] * XIn[iIn,jIn,kIn,:]
        end
      end
    end
  end  
  for kIn = 1 : nzIn
    for jIn = 1 : nyIn
      for jOut = 1 : nyOut
        for iOut = 1 : nxOut
          @views @.  Buf2[iOut,jOut,kIn,:] = Buf2[iOut,jOut,kIn,:] + 
            DG.IntXE2F[jOut,jIn] * Buf1[iOut,jIn,kIn,:]
        end
      end
    end
  end  
  for kIn = 1 : nzIn
    for kOut = 1 : nzOut
      for jOut = 1 : nyOut
        for iOut = 1 : nxOut
          @views @. XOut[iOut,jOut,kOut,:] = XOut[iOut,jOut,kOut,:] + 
            DG.IntZE2F[kOut,kIn] * Buf2[iOut,jOut,kIn,:]
        end
      end
    end
  end  
end
