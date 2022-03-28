using NLsolve

function SetIntersect(Circle1,Circle2,t)
function Intersect(t::Array{Float64,1})
  t1=t[1]
  t2=t[2]
  X1=t1*Circle1[:,1]+(1-t1)*Circle1[:,2]
  X2=t2*Circle2[:,1]+(1-t2)*Circle2[:,2]
  X1=X1/norm(X1)
  X2=X2/norm(X2)
  return [X1[1]-X2[1],X1[2]-X2[2]]
end  
end
function FindIntersectionPoint(Circle1,Circle2)
  t=[.5;.5]
  F=SetIntersect(Circle1,Circle2,t)
  res=nlsolve(F,t)
  return (res.zero)
end

function SetDistance(F,lon0,lat0,x)
function Distance(x::Array{Float64,1})
  ksi=x[1]
  eta=x[2]
  X=zeros(1,3);
  X[1:3]=0.25*((1-ksi)*(1-eta)*F.P[1:3,1]+
   (1+ksi)*(1-eta)*F.P[1:3,2]+
   (1+ksi)*(1+eta)*F.P[1:3,3]+
   (1-ksi)*(1+eta)*F.P[1:3,4]);
  r=norm(X);
  X=X/r;
  X0=sphereDeg2cart(lon0,lat0,1.0)
  return [X[1]-X0[1],X[2]-X0[2]]
end
end

function FindPointInCell(lon,lat,Grid)
  iFC = 0
  for iF=1:Grid.NumFaces
    @show iF  
    erg = PointInCell(lon,lat,Grid.Faces[iF])
    if erg
      iFC = iF  
      break
    end    
  end
  @show iFC
  x=zeros(2);
  F=SetDistance(Grid.Faces[iFC],lon,lat,x)
  res=nlsolve(F,x)
  @show res.zero
  println("AAA ",AAA)
  return (iFC,res.zero)
end

function PointInCell(lon,lat,Face)
  X=0.25*(Face.P[1:3,1]+ Face.P[1:3,2]+ Face.P[1:3,3]+ Face.P[1:3,4]);
  P=sphere2cart(lon,lat,1.0)
  t=FindIntersectionPoint([X P],[Face.P[1:3,1] Face.P[1:3,2]])
  @show t
  if 0<=t[1] && t[1]<=1 && 0<=t[2] && t[2]<=1
    return false  
  end  
  t=FindIntersectionPoint([X P],[Face.P[1:3,2] Face.P[1:3,3]])
  @show t
  if 0<=t[1] && t[1]<=1 && 0<=t[2] && t[2]<=1
    return false  
  end  
  t=FindIntersectionPoint([X P],[Face.P[1:3,3] Face.P[1:3,4]])
  @show t
  if 0<=t[1] && t[1]<=1 && 0<=t[2] && t[2]<=1
    return false  
  end  
  t=FindIntersectionPoint([X P],[Face.P[1:3,4] Face.P[1:3,1]])
  @show t
  if 0<=t[1] && t[1]<=1 && 0<=t[2] && t[2]<=1
    return false  
  end  
  return true
end

