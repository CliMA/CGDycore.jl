function triSolve!(x,tri,b)
# Thomas algorithm, b and x used as intermediate storage
# b is destroyed after the elemination
n = size(tri,2)
x[1] = tri[2,1]
@inbounds @fastmath for i=1:n-1
  x[i+1] = tri[2,i+1]-tri[3,i]/x[i]*tri[1,i+1]  
  b[i+1] = b[i+1]-tri[3,i]/x[i]*b[i]  
end
x[n] = b[n]/x[n]
@inbounds @fastmath for i=n-1:-1:1
  x[i]=(b[i] - tri[1,i+1]*x[i+1])/x[i]  
end
end

function mulUL!(tri,biU,biL)
n = size(biL,2)
tri[2,1] = tri[2,1] - biU[2,1]*biL[1,1] - biU[1,2]*biL[2,1]
tri[3,1] = tri[3,1] - biU[2,2]*biL[2,1]
@inbounds @fastmath for i=2:n-1
  tri[1,i] = tri[1,i] - biU[1,i]*biL[1,i]
  tri[2,i] = tri[2,i] - biU[2,i]*biL[1,i] - biU[1,i+1]*biL[2,i]
  tri[3,i] = tri[3,i] - biU[2,i+1]*biL[2,i]
end
i = n
tri[1,i] = tri[1,i] - biU[1,i]*biL[1,i]
tri[2,i] = tri[2,i] - biU[2,i]*biL[1,i]
end

function mulbiUv!(u,biU,v)
n = size(biU,2)
@inbounds @fastmath for i=1:n-1
  u[i] = u[i] + biU[2,i]*v[i] + biU[1,i+1]*v[i+1]
end
u[n] = u[n] + biU[2,n]*v[n]
end

function mulbiLv!(u,biL,v)
n = size(biL,2)
u[1] = u[1] + biL[1,1]*v[1]
@inbounds @fastmath for i=2:n
  u[i] = u[i] + biL[2,i-1]*v[i-1] + biL[1,i]*v[i]
end
end


function SchurSolve!(k,v,J,fac,Param)
n1=size(v,1);
n2=size(v,2);
n=n1*n2;
tri=J.tri
JWRho=J.JWRho
JRhoW=J.JRhoW
JWTh=J.JWTh
JThW=J.JThW
JWW=J.JWW
rRho=J.rRho
rTh=J.rTh
rw=J.rw
sw=J.sw

@views rRho=reshape(permute(v[:,:,1],[2,1]),n,1);
@views rTh=reshape(permute(v[:,:,5],[2,1]),n,1);
@views rw=reshape(permute(v[:,:,4],[2,1]),n,1);
invfac=1/fac;
invfac2=invfac/fac;
if Param.Damping
# sw=(spdiags(repmat(invfac2,n,1),0,n,n)-invfac*JWW-JWRho*JRhoW-JWTh*JThW)\
#   (invfac*rw+JWRho*rRho+JWTh*rTh);
  if Param.CompTri
    @views tri[1,:] .= 0
    @views tri[2,:] .= invfac2 .- invfac .* JWW[1,:]
    @views tri[3,:] .= 0
    mulUL!(tri,JWRho,JRhoW)
    mulUL!(tri,JWTh,JThW)
    Param.CompTri=false
  end
  rw .= invfac .* rw
  mulbiUv!(rw,JWRho,rRho)
  mulbiUv!(rw,JWTh,rTh)
  triSolve!(sw,tri,rw)
else
  if Param.CompTri  
    @views tri[1,:] .= 0
    @views tri[2,:] .= invfac2   
    @views tri[3,:] .= 0
    mulUL!(tri,JWRho,JRhoW)
    mulUL!(tri,JWTh,JThW)
    Param.CompTri=false
  end
  rw .= invfac .* rw
  mulbiUv!(rw,JWRho,rRho)
  mulbiUv!(rw,JWTh,rTh)
  triSolve!(sw,tri,rw)
end
mulbiLv!(rRho,JRhoW,sw)
mulbiLv!(rTh,JThW,sw)

@views k[:,:,1] .= fac .* permute(reshape(rRho,n2,n1),[2,1]);
@views k[:,:,2:3] .= fac*v[:,:,2:3];
@views k[:,:,4] .= permute(reshape(sw,n2,n1),[2,1]);
@views k[:,:,5] .= fac .* permute(reshape(rTh,n2,n1),[2,1]);
end


