function SetDistance(F,lon0,lat0,x)
function Distance(x::Array{Float64,1})
  ksi=x[1]
  eta=x[2]
  X=zeros(1,3);
  @show ksi
  @show eta
  X[1:3]=0.25*((1-ksi)*(1-eta)*F.P[1:3,1]+
   (1+ksi)*(1-eta)*F.P[1:3,2]+
   (1+ksi)*(1+eta)*F.P[1:3,3]+
   (1-ksi)*(1+eta)*F.P[1:3,4]);
  r=norm(X);
  X=X/r;
  X0=sphere2cart(lon0,lat0,1.0)
  @show X
  @show X0
  return [X[1]-X0[1],X[2]-X0[2]]
end
end

function FindPointInCell(lon,lat,Grid)
  iFC = 0
  for iF=1:Grid.NumFaces
    erg = PointInCell(lon,lat,Grid.Faces[iF])
    if erg
      iFC = iF  
      break
    end    
  end
  x=zeros(2);
  F=SetDistance(Grid.Faces[iFC],lon,lat,x)
  res=nlsolve(F,x)
  @show iFC
  @show res.zero
  println("AAA ",AAA)
  return (iFC,res.zero)
end

function PointInCell(lon,lat,Face)
  (lon1,lat1, _)=cart2sphere(Face.P[1,1],Face.P[2,1],Face.P[3,1]);
  (lon2,lat2, _)=cart2sphere(Face.P[1,2],Face.P[2,2],Face.P[3,2]);
  (lon3,lat3, _)=cart2sphere(Face.P[1,3],Face.P[2,3],Face.P[3,3]);
  (lon4,lat4, _)=cart2sphere(Face.P[1,4],Face.P[2,4],Face.P[3,4]);
  lonMin=min(lon1,lon2,lon3,lon4)
  if lon1-lonMin>pi
    lon1=2*pi -lon1
  end  
  if lon2-lonMin>pi
    lon2=2*pi -lon2
  end  
  if lon3-lonMin>pi
    lon3=2*pi -lon3
  end  
  if lon4-lonMin>pi
    lon4=2*pi -lon4
  end  
      
  poly = SVector{2,Float64}[(lon1,lat1),(lon2,lat2),(lon3,lat3),(lon4,lat4),(lon1,lat1)]
  algo=HaoSun()
  in=inpolygon(SVector(lon,lat),poly,algo)
  @show in
  if abs(in)>0
    @show (lon,lat)  
    @show (lon1,lat1)  
    @show (lon2,lat2)  
    @show (lon3,lat3)  
    @show (lon4,lat4)  
    return true
  else
    return false
  end  
end

