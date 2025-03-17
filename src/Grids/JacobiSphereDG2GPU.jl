function JacobiSphereDG2GPU!(X,dXdx,dXdxI,J,Rotate,FE,F,Rad)

  backend = get_backend(X)
  FT = eltype(X)

  NF = size(X,5)
  N = size(FE.xw,1)

  NFG = min(div(512,N*N*2),NF)
  group = (N, N, NFG)
  ndrange = (N, N, NF)

  KJacobiSphereDG2Kernel! = JacobiSphereDG2Kernel!(backend,group)

  KJacobiSphereDG2Kernel!(X,dXdx,dXdxI,J,Rotate,FE.xw,F,Rad,ndrange=ndrange)
end

@kernel inbounds = true function JacobiSphereDG2Kernel!(X,dXdx,dXdxI,JJ,Rotate,@Const(ksi),@Const(F),Rad)

  I, J, iF   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  FaceTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  NF = @uniform @ndrange()[3]

  dXdxLoc = @localmem eltype(X) (N,N,3,2,FaceTilesDim)

  eta = ksi
  if IF <= NF
    ID = I + (J - 1) * N
    ksi1 = ksi[I]
    ksi2 = ksi[J]
    X1 = eltype(X)(1/4) * (F[1,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,1,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,1,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    X2 = eltype(X)(1/4) * (F[1,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,2,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,2,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    X3 = eltype(X)(1/4) * (F[1,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)-ksi2) +
     F[2,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)-ksi2) +
     F[3,3,IF] * (eltype(X)(1)+ksi1)*(eltype(X)(1)+ksi2) +
     F[4,3,IF] * (eltype(X)(1)-ksi1)*(eltype(X)(1)+ksi2))
    r = sqrt(X1 * X1 + X2 * X2 + X3 * X3)

    f = Rad / r^3;
    dX1dX1 = f * (X2^2 + X3^2);
    dX1dX2 = -f * X1 * X2;
    dX1dX3 = -f * X1 * X3;
    dX2dX1 = dX1dX2;
    dX2dX2 = f *(X1^2 + X3^2);
    dX2dX3 = -f * X2 * X3;
    dX3dX1 = dX1dX3;
    dX3dX2 = dX2dX3;
    dX3dX3 = f * (X1^2 + X2^2);
    J1 = @SArray([dX1dX1    dX1dX2     dX1dX3
                  dX2dX1    dX2dX2     dX2dX3
                  dX3dX1    dX3dX2     dX3dX3])

    X1 = X1 / r
    X2 = X2 / r
    X3 = X3 / r

    B = @SArray([F[1,1,IF] F[2,1,IF] F[3,1,IF] F[4,1,IF]
                 F[1,2,IF] F[2,2,IF] F[3,2,IF] F[4,2,IF]
                 F[1,3,IF] F[2,3,IF] F[3,3,IF] F[4,3,IF]])
    C = @SArray([eltype(X)(-1)+ksi2  eltype(X)(-1)+ksi1
                 eltype(X)(1)-ksi2  eltype(X)(-1)-ksi1
                 eltype(X)(1)+ksi2   eltype(X)(1)+ksi1
                 eltype(X)(-1)-ksi2   eltype(X)(1)-ksi1])
    @views dXdxLoc[I,J,:,1:2,iF] .= eltype(X)(0.25) * J1 * B * C
    @views @. dXdx[:,1:2,1,ID,1,IF] = dXdxLoc[I,J,:,1:2,iF]

    @views JJ[ID,1,1,IF] = Determinant(dXdxLoc[I,J,:,1,iF],dXdxLoc[I,J,:,2,iF])
    @views pinvJac(dXdxI[:,:,1,ID,1,IF],dXdxLoc[I,J,:,:,iF]) 
    @views @. dXdxI[:,:,1,ID,1,IF] *= JJ[ID,1,1,IF]
    X[ID,1,1,1,IF] = X1 * Rad
    X[ID,1,2,1,IF] = X2 * Rad
    X[ID,1,3,1,IF] = X3 * Rad

    lon,lat,_ = cart2sphere(X1,X2,X3)
    Rotate[1,1,1,ID,1,IF] = -sin(lon)
    Rotate[2,1,1,ID,1,IF] = -sin(lat)*cos(lon)

    Rotate[1,2,1,ID,1,IF] = cos(lon)
    Rotate[2,2,1,ID,1,IF] = -sin(lat)*sin(lon)

    Rotate[1,3,1,ID,1,IF] = 0.0
    Rotate[2,3,1,ID,1,IF] = cos(lat)
  end  
end  

@inline function Determinant(a,b)

  D = (a[2]*b[3]-a[3]*b[2])^2 +
      (a[1]*b[3]-a[3]*b[1])^2 +
      (a[1]*b[2]-a[2]*b[1])^2
  D = sqrt(D)
end  

@inline function pinvJac(pinvD,D)
  g11 = D[1,1] * D[1,1] + D[2,1] * D[2,1] + D[3,1] * D[3,1]
  g12 = D[1,1] * D[1,2] + D[2,1] * D[2,2] + D[3,1] * D[3,2]
  g22 = D[1,2] * D[1,2] + D[2,2] * D[2,2] + D[3,2] * D[3,2]
  det = g11 * g22 - g12^2
  i11 = g22 / det
  i21 = -g12 / det
  i12 = -g12 / det
  i22 = g11 / det
  pinvD[1,1] = D[1,1] * i11 + D[1,2] * i21
  pinvD[2,1] = D[1,1] * i12 + D[1,2] * i22
  pinvD[1,2] = D[2,1] * i11 + D[2,2] * i21
  pinvD[2,2] = D[2,1] * i12 + D[2,2] * i22
  pinvD[1,3] = D[3,1] * i11 + D[3,2] * i21
  pinvD[2,3] = D[3,1] * i12 + D[3,2] * i22

end


