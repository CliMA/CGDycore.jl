@kernel inbounds = true function MomentumVectorInvariantCoriolisKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(Glob),CoriolisFun,GridForm,Curvature)

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  uCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+2)
  tempuZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempvZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempwZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  YC = @localmem eltype(F) (N,N,3,2,2,ColumnTilesDim)

  if Iz <= Nz
    uCol[I,J,iz] = U[1,Iz,ind,2]
    vCol[I,J,iz] = U[1,Iz,ind,3]
    RhoCol[I,J,iz+1] = U[1,Iz,ind,1]
    wCol[I,J,iz+1] = U[1,Iz,ind,4]
    if iz == 1 && Iz == 1
      wCol[I,J,1] = -(dXdxI[3,1,1,ID,1,IF] * U[1,Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[1,Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      RhoCol[I,J,1] = U[1,Iz,ind,1]
    elseif iz == 1
      wCol[I,J,1] = U[1,Iz-1,ind,4]
      RhoCol[I,J,1] = U[1,Iz-1,ind,1]
    elseif iz == ColumnTilesDim && Iz < Nz  
      RhoCol[I,J,ColumnTilesDim+2] = U[1,Iz+1,ind,1]
    end   
  end 
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    # Dz*(dx33*v - dx32*w)
    tempuZ[I,J,1,iz+1] = dXdxI[3,3,1,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[3,2,1,ID,Iz,IF] * wCol[I,J,iz]
    tempuZ[I,J,2,iz+1] = dXdxI[3,3,2,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[3,2,2,ID,Iz,IF] * wCol[I,J,iz+1]
    # Dz*(dx33*u - dx31*w)
    tempvZ[I,J,1,iz+1] = dXdxI[3,3,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,1,ID,Iz,IF] * wCol[I,J,iz]
    tempvZ[I,J,2,iz+1] = dXdxI[3,3,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,2,ID,Iz,IF] * wCol[I,J,iz+1]
    # Dz*(dx32*u - dx31*v)
    tempwZ[I,J,1,iz+1] = dXdxI[3,2,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,1,ID,Iz,IF] * vCol[I,J,iz]
    tempwZ[I,J,2,iz+1] = dXdxI[3,2,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,2,ID,Iz,IF] * vCol[I,J,iz] 
    if iz == 1 && Iz > 1
      um1 = U[1,Iz-1,ind,2]  
      vm1 = U[1,Iz-1,ind,3]  
      tempuZ[I,J,2,1] = dXdxI[3,3,2,ID,Iz-1,IF] * vm1 - dXdxI[3,2,2,ID,Iz-1,IF] * wCol[I,J,1]
      tempvZ[I,J,2,1] = dXdxI[3,3,2,ID,Iz-1,IF] * um1 - dXdxI[3,1,2,ID,Iz-1,IF] * wCol[I,J,1]
      tempwZ[I,J,2,1] = dXdxI[3,2,2,ID,Iz-1,IF] * um1 - dXdxI[3,1,2,ID,Iz-1,IF] * vm1 
    end  
    if iz == ColumnTilesDim && Iz < Nz
      up1 = U[1,Iz+1,ind,2]  
      vp1 = U[1,Iz+1,ind,3]  
      tempuZ[I,J,1,iz+2] = dXdxI[3,3,1,ID,Iz+1,IF] * vp1 - dXdxI[3,2,1,ID,Iz+1,IF] * wCol[I,J,iz+1]
      tempvZ[I,J,1,iz+2] = dXdxI[3,3,1,ID,Iz+1,IF] * up1 - dXdxI[3,1,1,ID,Iz+1,IF] * wCol[I,J,iz+1]
      tempwZ[I,J,1,iz+2] = dXdxI[3,2,1,ID,Iz+1,IF] * up1 - dXdxI[3,1,1,ID,Iz+1,IF] * vp1 
    end  
  end  
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz

# uDot = - v*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))  
#        - w*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w))
# vDot =   u*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))
#        - w*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
# wDot =   u*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)) 
#          v*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
#   U = Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w)    
#   V = Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)
#   W = Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v)
#   (U,V,W) = (-y1,y2,-y3) with y = dx_a x (u,v,w). The horizontal derivatives act on R^T y,
#   the Cartesian components, and the result is rotated back with R (rows e_lam, e_phi, e_r).
#   This is the curl of the Cartesian velocity, R is constant along the vertical.
    R = CurlFrame(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],GridForm)
    @unroll for K = 1 : 2
      wK = wCol[I,J,iz+K-1]
      @unroll for a = 1 : 2
        y1 = dXdxI[a,2,K,ID,Iz,IF] * wK - dXdxI[a,3,K,ID,Iz,IF] * vCol[I,J,iz]
        y2 = dXdxI[a,3,K,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[a,1,K,ID,Iz,IF] * wK
        y3 = dXdxI[a,1,K,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[a,2,K,ID,Iz,IF] * uCol[I,J,iz]
        YC[I,J,1,K,a,iz] = R[1,1] * y1 + R[2,1] * y2 + R[3,1] * y3
        YC[I,J,2,K,a,iz] = R[1,2] * y1 + R[2,2] * y2 + R[3,2] * y3
        YC[I,J,3,K,a,iz] = R[1,3] * y1 + R[2,3] * y2 + R[3,3] * y3
      end
    end
    
  end
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    FU = eltype(F)(0)  
    FV = eltype(F)(0)  
    FW1 = eltype(F)(0)  
    FW2 = eltype(F)(0)  
    if Iz > 1
      FluxUZ = eltype(F)(0.5) * (tempuZ[I,J,1,iz+1] - tempuZ[I,J,2,iz])
      FluxVZ = eltype(F)(0.5) * (tempvZ[I,J,1,iz+1] - tempvZ[I,J,2,iz])
      FluxWZ = eltype(F)(0.5) * (tempwZ[I,J,1,iz+1] - tempwZ[I,J,2,iz])
      FU += (-vCol[I,J,iz] * FluxWZ - wCol[I,J,iz] * FluxVZ)
      FV += (uCol[I,J,iz] * FluxWZ - wCol[I,J,iz] * FluxUZ)
      FW1 += ( uCol[I,J,iz] * FluxVZ + vCol[I,J,iz] * FluxUZ)
    end
    if Iz < Nz
      FluxUZ = eltype(F)(0.5) * (tempuZ[I,J,1,iz+2] - tempuZ[I,J,2,iz+1])
      FluxVZ = eltype(F)(0.5) * (tempvZ[I,J,1,iz+2] - tempvZ[I,J,2,iz+1])
      FluxWZ = eltype(F)(0.5) * (tempwZ[I,J,1,iz+2] - tempwZ[I,J,2,iz+1])
      FU += (-vCol[I,J,iz] * FluxWZ - wCol[I,J,iz+1] * FluxVZ)
      FV += ( uCol[I,J,iz] * FluxWZ - wCol[I,J,iz+1] * FluxUZ)
      FW2 += ( uCol[I,J,iz] * FluxVZ + vCol[I,J,iz] * FluxUZ)
    end
    U1 = eltype(F)(0.5) * (tempuZ[I,J,2,iz+1] - tempuZ[I,J,1,iz+1])
    U2 = eltype(F)(0.5) * (tempuZ[I,J,2,iz+1] - tempuZ[I,J,1,iz+1])
    V1 = eltype(F)(0.5) * (tempvZ[I,J,2,iz+1] - tempvZ[I,J,1,iz+1])
    V2 = eltype(F)(0.5) * (tempvZ[I,J,2,iz+1] - tempvZ[I,J,1,iz+1])
    W1 = eltype(F)(0.5) * (tempwZ[I,J,2,iz+1] - tempwZ[I,J,1,iz+1])
    W2 = eltype(F)(0.5) * (tempwZ[I,J,2,iz+1] - tempwZ[I,J,1,iz+1])
    c11 = D[I,1] * YC[1,J,1,1,1,iz] + D[J,1] * YC[I,1,1,1,2,iz]
    c21 = D[I,1] * YC[1,J,2,1,1,iz] + D[J,1] * YC[I,1,2,1,2,iz]
    c31 = D[I,1] * YC[1,J,3,1,1,iz] + D[J,1] * YC[I,1,3,1,2,iz]
    c12 = D[I,1] * YC[1,J,1,2,1,iz] + D[J,1] * YC[I,1,1,2,2,iz]
    c22 = D[I,1] * YC[1,J,2,2,1,iz] + D[J,1] * YC[I,1,2,2,2,iz]
    c32 = D[I,1] * YC[1,J,3,2,1,iz] + D[J,1] * YC[I,1,3,2,2,iz]
    for k = 2 : N
      c11 += D[I,k] * YC[k,J,1,1,1,iz] + D[J,k] * YC[I,k,1,1,2,iz]
      c21 += D[I,k] * YC[k,J,2,1,1,iz] + D[J,k] * YC[I,k,2,1,2,iz]
      c31 += D[I,k] * YC[k,J,3,1,1,iz] + D[J,k] * YC[I,k,3,1,2,iz]
      c12 += D[I,k] * YC[k,J,1,2,1,iz] + D[J,k] * YC[I,k,1,2,2,iz]
      c22 += D[I,k] * YC[k,J,2,2,1,iz] + D[J,k] * YC[I,k,2,2,2,iz]
      c32 += D[I,k] * YC[k,J,3,2,1,iz] + D[J,k] * YC[I,k,3,2,2,iz]
    end  
    R = CurlFrame(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],GridForm)
    U1 -= R[1,1] * c11 + R[1,2] * c21 + R[1,3] * c31
    U2 -= R[1,1] * c12 + R[1,2] * c22 + R[1,3] * c32
    V1 += R[2,1] * c11 + R[2,2] * c21 + R[2,3] * c31
    V2 += R[2,1] * c12 + R[2,2] * c22 + R[2,3] * c32
    W1 -= R[3,1] * c11 + R[3,2] * c21 + R[3,3] * c31
    W2 -= R[3,1] * c12 + R[3,2] * c22 + R[3,3] * c32

#   Shallow atmosphere: with the shallow metric the curl of the Cartesian velocity adds the
#   tendencies u*w/a, v*w/a, -(u^2+v^2)/a, which are subtracted again.
#   Energy neutral in each half cell.
    w1 = wCol[I,J,iz]
    w2 = wCol[I,J,iz+1]
    KinH = uCol[I,J,iz] * uCol[I,J,iz] + vCol[I,J,iz] * vCol[I,J,iz]
    FU -= Curvature * uCol[I,J,iz] * (JJ[ID,1,Iz,IF] * w1 + JJ[ID,2,Iz,IF] * w2)
    FV -= Curvature * vCol[I,J,iz] * (JJ[ID,1,Iz,IF] * w1 + JJ[ID,2,Iz,IF] * w2)
    FW1 += Curvature * JJ[ID,1,Iz,IF] * KinH
    FW2 += Curvature * JJ[ID,2,Iz,IF] * KinH

#   Coriolis, in each half cell with its vertical velocity (energy neutral also with
#   the non-traditional terms)
    FuCor, FvCor, FwCor = CoriolisFun(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],
      uCol[I,J,iz],vCol[I,J,iz],wCol[I,J,iz],wCol[I,J,iz])
    FU += FuCor * JJ[ID,1,Iz,IF]
    FV += FvCor * JJ[ID,1,Iz,IF]
    FW1 += FwCor * JJ[ID,1,Iz,IF]
    FuCor, FvCor, FwCor = CoriolisFun(X[ID,2,1,Iz,IF],X[ID,2,2,Iz,IF],X[ID,2,3,Iz,IF],
      uCol[I,J,iz],vCol[I,J,iz],wCol[I,J,iz+1],wCol[I,J,iz+1])
    FU += FuCor * JJ[ID,2,Iz,IF]
    FV += FvCor * JJ[ID,2,Iz,IF]
    FW2 += FwCor * JJ[ID,2,Iz,IF]

    @atomic :monotonic F[1,Iz,ind,2] += (-vCol[I,J,iz] * W1 - wCol[I,J,iz] * V1 -
      vCol[I,J,iz] * W2 - wCol[I,J,iz+1] * V2 + FU) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[1,Iz,ind,3] += (uCol[I,J,iz] * W1 - wCol[I,J,iz] * U1 +
      uCol[I,J,iz] * W2 - wCol[I,J,iz+1] * U2 + FV) / (M[Iz,ind,1] + M[Iz,ind,2])
    if Iz > 1  
      @atomic :monotonic F[1,Iz-1,ind,4] += RhoCol[I,J,iz+1] * (uCol[I,J,iz] * V1 + vCol[I,J,iz] * U1 + FW1) /
        (RhoCol[I,J,iz] * M[Iz-1,ind,2] + RhoCol[I,J,iz+1] * M[Iz,ind,1])
    end
    if Iz < Nz
      @atomic :monotonic F[1,Iz,ind,4] += RhoCol[I,J,iz+1] * (uCol[I,J,iz] * V2 + vCol[I,J,iz] * U2 + FW2) / 
        (RhoCol[I,J,iz+2] * M[Iz+1,ind,1] + RhoCol[I,J,iz+1] * M[Iz,ind,2])
    end  
  end
end

# Rotation to the spherical frame at a node, rows e_lam, e_phi, e_r
@inline function CurlFrame(x,y,z,::Grids.SphericalGrid)
  lon,lat,_ = Grids.cart2sphere(x,y,z)
  return Grids.MCart2Sphere(lon,lat)
end

@inline function CurlFrame(x,y,z,::Grids.CartesianGrid)
  FT = eltype(x)
  return @SMatrix [FT(1) FT(0) FT(0); FT(0) FT(1) FT(0); FT(0) FT(0) FT(1)]
end

# Coefficient of the curvature terms that are removed for the shallow atmosphere
CurvatureCoefficient(Rad,::Grids.SphericalGrid,::Models.CompressibleDeep) = zero(Rad)
CurvatureCoefficient(Rad,::Grids.SphericalGrid,Equation) = one(Rad) / Rad
CurvatureCoefficient(Rad,::Grids.CartesianGrid,Equation) = zero(Rad)

@kernel inbounds = true function MomentumVectorInvariantCoriolisDraftKernel!(F,@Const(U),@Const(w),@Const(Rho),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(Glob),CoriolisFun,GridForm,Curvature)

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  uCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  tempuZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempvZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  tempwZ = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  YC = @localmem eltype(F) (N,N,3,2,2,ColumnTilesDim)

  if Iz <= Nz
    uCol[I,J,iz] = U[1,Iz,ind,2]
    vCol[I,J,iz] = U[1,Iz,ind,3]
    wCol[I,J,iz+1] = w[Iz,ind,IE]
    RhoCol[I,J,iz+1] = Rho[Iz,ind,IE]
    if iz == 1 && Iz == 1
      wCol[I,J,1] = -(dXdxI[3,1,1,ID,1,IF] * U[1,Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[1,Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      RhoCol[I,J,1] = Rho[1,ind,IE]
    elseif iz == 1
      wCol[I,J,1] = w[Iz-1,ind,IE]
      RhoCol[I,J,1] = Rho[Iz-1,ind,IE]
    end   
  end 
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    # Dz*(dx33*v - dx32*w)
    tempuZ[I,J,1,iz+1] = dXdxI[3,3,1,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[3,2,1,ID,Iz,IF] * wCol[I,J,iz]
    tempuZ[I,J,2,iz+1] = dXdxI[3,3,2,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[3,2,2,ID,Iz,IF] * wCol[I,J,iz+1]
    # Dz*(dx33*u - dx31*w)
    tempvZ[I,J,1,iz+1] = dXdxI[3,3,1,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,1,ID,Iz,IF] * wCol[I,J,iz]
    tempvZ[I,J,2,iz+1] = dXdxI[3,3,2,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[3,1,2,ID,Iz,IF] * wCol[I,J,iz+1]
    if iz == 1 && Iz > 1
      um1 = U[1,Iz-1,ind,2]  
      vm1 = U[1,Iz-1,ind,3]  
      tempuZ[I,J,2,1] = dXdxI[3,3,2,ID,Iz-1,IF] * vm1 - dXdxI[3,2,2,ID,Iz-1,IF] * wCol[I,J,1]
      tempvZ[I,J,2,1] = dXdxI[3,3,2,ID,Iz-1,IF] * um1 - dXdxI[3,1,2,ID,Iz-1,IF] * wCol[I,J,1]
    end  
    if iz == ColumnTilesDim && Iz < Nz
      up1 = U[1,Iz+1,ind,2]  
      vp1 = U[1,Iz+1,ind,3]  
      tempuZ[I,J,1,iz+2] = dXdxI[3,3,1,ID,Iz+1,IF] * vp1 - dXdxI[3,2,1,ID,Iz+1,IF] * wCol[I,J,iz+1]
      tempvZ[I,J,1,iz+2] = dXdxI[3,3,1,ID,Iz+1,IF] * up1 - dXdxI[3,1,1,ID,Iz+1,IF] * wCol[I,J,iz+1]
    end  
  end  
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz

# uDot = - v*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))  
#        - w*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w))
# vDot =   u*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))
#        - w*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
# wDot =   u*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)) 
#          v*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
#   U = Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w)    
#   V = Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)
#   W = Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v)
#   Horizontal derivatives of the Cartesian components, see MomentumVectorInvariantCoriolisKernel!
    R = CurlFrame(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],GridForm)
    @unroll for K = 1 : 2
      wK = wCol[I,J,iz+K-1]
      @unroll for a = 1 : 2
        y1 = dXdxI[a,2,K,ID,Iz,IF] * wK - dXdxI[a,3,K,ID,Iz,IF] * vCol[I,J,iz]
        y2 = dXdxI[a,3,K,ID,Iz,IF] * uCol[I,J,iz] - dXdxI[a,1,K,ID,Iz,IF] * wK
        y3 = dXdxI[a,1,K,ID,Iz,IF] * vCol[I,J,iz] - dXdxI[a,2,K,ID,Iz,IF] * uCol[I,J,iz]
        YC[I,J,1,K,a,iz] = R[1,1] * y1 + R[2,1] * y2 + R[3,1] * y3
        YC[I,J,2,K,a,iz] = R[1,2] * y1 + R[2,2] * y2 + R[3,2] * y3
        YC[I,J,3,K,a,iz] = R[1,3] * y1 + R[2,3] * y2 + R[3,3] * y3
      end
    end
    
  end
  @synchronize 

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    FW1 = eltype(F)(0)  
    FW2 = eltype(F)(0)  
    if Iz > 1
      FluxUZ = eltype(F)(0.5) * (tempuZ[I,J,1,iz+1] - tempuZ[I,J,2,iz])
      FluxVZ = eltype(F)(0.5) * (tempvZ[I,J,1,iz+1] - tempvZ[I,J,2,iz])
      FW1 += ( uCol[I,J,iz] * FluxVZ + vCol[I,J,iz] * FluxUZ)
    end
    if Iz < Nz
      FluxUZ = eltype(F)(0.5) * (tempuZ[I,J,1,iz+2] - tempuZ[I,J,2,iz+1])
      FluxVZ = eltype(F)(0.5) * (tempvZ[I,J,1,iz+2] - tempvZ[I,J,2,iz+1])
      FW2 += ( uCol[I,J,iz] * FluxVZ + vCol[I,J,iz] * FluxUZ)
    end
    U1 = eltype(F)(0.5) * (tempuZ[I,J,2,iz+1] - tempuZ[I,J,1,iz+1])
    U2 = eltype(F)(0.5) * (tempuZ[I,J,2,iz+1] - tempuZ[I,J,1,iz+1])
    V1 = eltype(F)(0.5) * (tempvZ[I,J,2,iz+1] - tempvZ[I,J,1,iz+1])
    V2 = eltype(F)(0.5) * (tempvZ[I,J,2,iz+1] - tempvZ[I,J,1,iz+1])
    c11 = D[I,1] * YC[1,J,1,1,1,iz] + D[J,1] * YC[I,1,1,1,2,iz]
    c21 = D[I,1] * YC[1,J,2,1,1,iz] + D[J,1] * YC[I,1,2,1,2,iz]
    c31 = D[I,1] * YC[1,J,3,1,1,iz] + D[J,1] * YC[I,1,3,1,2,iz]
    c12 = D[I,1] * YC[1,J,1,2,1,iz] + D[J,1] * YC[I,1,1,2,2,iz]
    c22 = D[I,1] * YC[1,J,2,2,1,iz] + D[J,1] * YC[I,1,2,2,2,iz]
    c32 = D[I,1] * YC[1,J,3,2,1,iz] + D[J,1] * YC[I,1,3,2,2,iz]
    for k = 2 : N
      c11 += D[I,k] * YC[k,J,1,1,1,iz] + D[J,k] * YC[I,k,1,1,2,iz]
      c21 += D[I,k] * YC[k,J,2,1,1,iz] + D[J,k] * YC[I,k,2,1,2,iz]
      c31 += D[I,k] * YC[k,J,3,1,1,iz] + D[J,k] * YC[I,k,3,1,2,iz]
      c12 += D[I,k] * YC[k,J,1,2,1,iz] + D[J,k] * YC[I,k,1,2,2,iz]
      c22 += D[I,k] * YC[k,J,2,2,1,iz] + D[J,k] * YC[I,k,2,2,2,iz]
      c32 += D[I,k] * YC[k,J,3,2,1,iz] + D[J,k] * YC[I,k,3,2,2,iz]
    end  
    R = CurlFrame(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],GridForm)
    U1 -= R[1,1] * c11 + R[1,2] * c21 + R[1,3] * c31
    U2 -= R[1,1] * c12 + R[1,2] * c22 + R[1,3] * c32
    V1 += R[2,1] * c11 + R[2,2] * c21 + R[2,3] * c31
    V2 += R[2,1] * c12 + R[2,2] * c22 + R[2,3] * c32

#   Shallow atmosphere: subtract the tendency -(u^2+v^2)/a, see MomentumVectorInvariantCoriolisKernel!
    KinH = uCol[I,J,iz] * uCol[I,J,iz] + vCol[I,J,iz] * vCol[I,J,iz]
    FW1 += Curvature * JJ[ID,1,Iz,IF] * KinH
    FW2 += Curvature * JJ[ID,2,Iz,IF] * KinH

#   Coriolis, in each half cell
    _, _, FwCor = CoriolisFun(X[ID,1,1,Iz,IF],X[ID,1,2,Iz,IF],X[ID,1,3,Iz,IF],
      uCol[I,J,iz],vCol[I,J,iz],wCol[I,J,iz],wCol[I,J,iz])
    FW1 += FwCor * JJ[ID,1,Iz,IF]
    _, _, FwCor = CoriolisFun(X[ID,2,1,Iz,IF],X[ID,2,2,Iz,IF],X[ID,2,3,Iz,IF],
      uCol[I,J,iz],vCol[I,J,iz],wCol[I,J,iz+1],wCol[I,J,iz+1])
    FW2 += FwCor * JJ[ID,2,Iz,IF]

    if Iz > 1  
      @atomic :monotonic F[1,Iz-1,ind,4] += RhoCol[I,J,iz] * (uCol[I,J,iz] * V1 + vCol[I,J,iz] * U1 + FW1) /
        (RhoCol[I,J,iz-1] * M[Iz-1,ind,2] + RhoCol[I,J,iz] * M[Iz,ind,1])
    end
    if Iz < Nz
      @atomic :monotonic F[1,Iz,ind,4] += RhoCol[I,J,iz] * (uCol[I,J,iz] * V2 + vCol[I,J,iz] * U2 + FW2) / 
        (RhoCol[I,J,iz+1] * M[Iz+1,ind,1] + RhoCol[I,J,iz] * M[Iz,ind,2])
    end  
  end
end

@kernel inbounds = true function VerticalDiffusionMomentumKernel!(F,@Const(U),@Const(K),
  @Const(dz))
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if Iz < Nz && IC <= NumG
    fac = eltype(F)(2) * K[Iz,IC] / (dz[Iz,IC] + dz[Iz+1,IC]) 
    facDiv1 = eltype(dz)(1) / dz[Iz,IC] / U[1,Iz,IC,1]
    facDiv2 = -eltype(dz)(1) / dz[Iz+1,IC] / U[1,Iz+1,IC,1]
    grad = fac * (U[1,Iz+1,IC,2] - U[1,Iz,IC,2])
    @atomic :monotonic F[1,Iz,IC,2] +=  facDiv1 * grad
    @atomic :monotonic F[1,Iz+1,IC,2] += facDiv2 * grad 
    grad = fac * (U[1,Iz+1,IC,3] - U[1,Iz,IC,3])
    @atomic :monotonic F[1,Iz,IC,3] +=  facDiv1 * grad
    @atomic :monotonic F[1,Iz+1,IC,3] += facDiv2 * grad 
  end  
end  


@kernel inbounds = true function SurfaceFluxMomentumKernel!(F,@Const(U),@Const(nS),@Const(CM),@Const(dz))
  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    uCol = U[1,1,IC,2]
    vCol = U[1,1,IC,3]
    WS = -(nS[1,IC]* uCol + nS[2,IC] * vCol) / nS[3,IC]
    ww = 0.5 * (WS + U[1,1,IC,4])
    nU = nS[1,IC] * uCol + nS[2,IC] * vCol + nS[3,IC] * ww
    uStar = sqrt((uCol - nS[1,IC] * nU)^2 + (vCol - nS[2,IC] * nU)^2 + (ww - nS[3,IC] * nU)^2)
    facDiv = -eltype(dz)(2) / dz[1,IC] * CM[IC] * uStar
    @atomic :monotonic F[1,1,IC,2] += facDiv * (uCol - nU * nS[1,IC])
    @atomic :monotonic F[1,1,IC,3] += facDiv * (vCol - nU * nS[2,IC])
  end
end  
