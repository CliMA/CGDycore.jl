function JacobiDG3(CG,F,z,Topo,Topography,zs)
ksi=CG.xe
eta=CG.xe
zeta=CG.xwZ
z1=z[1]
z2=z[2]
n=CG.OrdPoly+1
nz=CG.OrdPolyZ+1
XE=zeros(n,n,nz,3)
X=zeros(n,n,nz,3)
dXdx=zeros(n,n,nz,3,3)
dXdxI=zeros(3,3,nz,n,n)
J=zeros(n,n,nz)
theta=zeros(n,n)
(_,DS) = DG.DerivativeMatrixSingle(CG.OrdPoly)
(_,DSZ) = DG.DerivativeMatrixSingle(CG.OrdPolyZ)
@inbounds for j=1:n
  @inbounds for i=1:n
    @inbounds for k=1:nz
    XE[i,j,k,1]=0.25*((1-ksi[i])*(1-eta[j])*F.P[1].x+
                    (1+ksi[i])*(1-eta[j])*F.P[2].x+
                    (1+ksi[i])*(1+eta[j])*F.P[3].x+
                    (1-ksi[i])*(1+eta[j])*F.P[4].x)
    XE[i,j,k,2]=0.25*((1-ksi[i])*(1-eta[j])*F.P[1].y+
                    (1+ksi[i])*(1-eta[j])*F.P[2].y+
                    (1+ksi[i])*(1+eta[j])*F.P[3].y+
                    (1-ksi[i])*(1+eta[j])*F.P[4].y)
    z=0.5*((1-zeta[k])*z1+(1+zeta[k])*z2)
    (XE[i,j,k,3],dXdx[i,j,k,3,3])=Topo(XE[i,j,k,1],XE[i,j,k,2],z,z,Topography,zs)
    dXdx[i,j,k,3,3]=0.5*dXdx[i,j,k,3,3]*(z2-z1)
    dXdx[i,j,k,3,3]=0.0
    end
  end
end

  ChangeBasis!(X,XE,CG)
  @inbounds for k = 1 : nz
    @inbounds for j = 1 : n
      @inbounds for i = 1 : n
        @inbounds for l = 1 : n
          @views @. dXdx[i,j,k,:,1] = dXdx[i,j,k,:,1] + DS[i,l] * X[l,j,k,:]
        end
      end
    end
  end

  @inbounds for k = 1 : nz
    @inbounds for i = 1 : n
      @inbounds for j = 1 : n
        @inbounds for l = 1 : n
          @views @. dXdx[i,j,k,:,2] = dXdx[i,j,k,:,2] + DS[j,l] * X[i,l,k,:]
        end
      end
    end
  end
  @inbounds for j = 1 : n
    @inbounds for i = 1 : n
      @inbounds for k = 1 : nz
        @inbounds for l = 1 : nz
          @views @. dXdx[i,j,k,:,3] = dXdx[i,j,k,:,3] + DSZ[k,l] * X[i,j,l,:]
        end
      end
    end
  end

@inbounds for j=1:n
  @inbounds for i=1:n
    @inbounds for k=1:nz
      J[i,j,k]=det(reshape(dXdx[i,j,k,:,:],3,3))
      dXdxI[:,:,k,i,j]=inv(reshape(dXdx[i,j,k,:,:],3,3))*J[i,j,k]
    end
  end
end
X = reshape(X,n*n,nz,3)
J = reshape(J,n*n,nz)
dXdx = reshape(dXdx,n*n,nz,3,3)
dXdxI = reshape(dXdxI,3,3,nz,n*n)
return X,J,dXdx,dXdxI
end


function ChangeBasis!(XOut,XIn,CG)

  nxOut = size(XOut,1)
  nyOut = size(XOut,2)
  nzOut = size(XOut,3)
  nxIn = size(XIn,1)
  nyIn = size(XIn,2)
  nzIn = size(XIn,3)

  Buf1 = zeros(nxOut,nyIn,nzIn,3)
  Buf2 = zeros(nxOut,nyOut,nzIn,3)

  @inbounds for kIn = 1 : nzIn
    @inbounds for jIn = 1 : nyIn
      @inbounds for iIn = 1 : nxIn
        @inbounds for iOut = 1 : nxOut
          @views @.  Buf1[iOut,jIn,kIn,:] = Buf1[iOut,jIn,kIn,:] + 
            CG.IntXE2F[iOut,iIn] * XIn[iIn,jIn,kIn,:]
        end
      end
    end
  end  
  @inbounds for kIn = 1 : nzIn
    @inbounds for jIn = 1 : nyIn
      @inbounds for jOut = 1 : nyOut
        @inbounds for iOut = 1 : nxOut
          @views @.  Buf2[iOut,jOut,kIn,:] = Buf2[iOut,jOut,kIn,:] + 
            CG.IntXE2F[jOut,jIn] * Buf1[iOut,jIn,kIn,:]
        end
      end
    end
  end  
  @inbounds for kIn = 1 : nzIn
    @inbounds for kOut = 1 : nzOut
      @inbounds for jOut = 1 : nyOut
        @inbounds for iOut = 1 : nxOut
          @views @. XOut[iOut,jOut,kOut,:] = XOut[iOut,jOut,kOut,:] + 
            CG.IntZE2F[kOut,kIn] * Buf2[iOut,jOut,kIn,:]
        end
      end
    end
  end  
end
