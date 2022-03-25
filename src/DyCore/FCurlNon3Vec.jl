function FCurlNon3Vec!(FuHat,v1CG,v2CG,wCG,wCCG,CG,Param)
OP=CG.OrdPoly+1
NF=Param.Grid.NumFaces
nz=Param.Grid.nz
uPos=Param.uPos
vPos=Param.vPos
wPos = Param.wPos
JC = Param.cache.JC 

dXdxIC11 = Param.dXdxIC11
dXdxIC21 = Param.dXdxIC21
dXdxIC31 = Param.dXdxIC31
dXdxIC12 = Param.dXdxIC12
dXdxIC22 = Param.dXdxIC22
dXdxIC32 = Param.dXdxIC32
dXdxIC33 = Param.dXdxIC33

dXdxIC = Param.cache.dXdxIC
dXdxIF = Param.cache.dXdxIF

vHat = Param.CacheC1
Vort3 = Param.CacheC2
Temp = Param.CacheC3
DZvHat3 = Param.CacheC4

wHat = Param.CacheF1
DXwHat11 = Param.CacheF2
DYwHat21 = Param.CacheF3
DZuuHat31 = Param.CacheF4
DZwwHat31 = Param.CacheC5

DXwHat12 = Param.CacheF2
DYwHat22 = Param.CacheF3
DZuuHat32 = Param.CacheF4
DZwwHat32 = Param.CacheC5


# Fu(1)
# -vv\left(DX(a^{12}u-a^{11}v)+(a^{22}u-a^{21}v)DY^T
# +\frac{d(a^{32}u-a^{31}v)}{dz}\right)
#
# ww(DXa^{11}w+a^{21}wDY^T+\frac{da^{31}w}{dz}-\frac{da^{33}u}{dz})

# Fu(2)
# uu\left(DX(a^{12}u-a^{11}v)+(a^{22}u-a^{21}v)DY^T
# +\frac{d(a^{32}u-a^{31}v)}{dz}\right)
#
# ww(DXa^{12}w+a^{22}wDY^T+\frac{da^{32}w}{dz}-\frac{d a^{33}v}{dz})

# Fw
# -uu((DXa^{11}w+a^{21}wDY^T+\frac{da^{31}w}{dz}-\frac{da^{33}u}{dz})
# -vv*((DXa^{12}w+a^{22}wDY^T+\frac{da^{32}w}{dz}-\frac{da^{33}v}{dz})

vHat = dXdxIC12.*v1CG .- dXdxIC11.*v2CG
mul!(reshape(Vort3,OP,OP*NF*nz),CG.DS,reshape(vHat,OP,OP*nz*NF))
vHat = dXdxIC22.*v1CG .- dXdxIC21.*v2CG
mul!(reshape(PermutedDimsArray(Temp,(2,1,3,4)),OP,OP*NF*nz),CG.DS,reshape(PermutedDimsArray(vHat,(2,1,3,4)),OP,OP*nz*NF))
Vort3 .= Vort3 .+ Temp
vHat = dXdxIC32.*v1CG .- dXdxIC31.*v2CG
if nz>1
  @views Vort3[:,:,:,1] .+= 0.5 .* (vHat[:,:,:,2] .- vHat[:,:,:,1])
  @views Vort3[:,:,:,2:nz-1] .+= 0.25 .* (vHat[:,:,:,3:nz] .- vHat[:,:,:,1:nz-2])
  @views Vort3[:,:,:,nz] .+= 0.5 .* (vHat[:,:,:,nz] .- vHat[:,:,:,nz-1])
end
if Param.Coriolis
  str = Param.CoriolisType
  if str == "Sphere"
      lat = Param.cache.lat
      JC = Param.cache.JC
      Omega = Param.Omega
      @views Vort3 .= Vort3 .- reshape((2.0*Omega).* sin.(
        repmat(reshape(lat[:,:,:],OP*OP*NF,1),1,nz))
        ,OP,OP,NF,nz).* JC
  elseif str == "Beta-Plane"
      JC = Param.cache.JC
      X = Param.cache.X
      beta0 = Param.beta0
      f0 = Param.f0
      y0 = Param.y0
      @views Vort3 .= Vort3 .- (f0 .+ beta0 .* (
        reshape(abs.(0.5.*(X[:,:,1,2,:,:] .+ X[:,:,2,2,:,:])),OP,OP,NF,nz) .- y0)) .* JC
  end
end
@views FuHat[:,:,:,:,uPos] .-= Vort3.*v2CG 
@views FuHat[:,:,:,:,vPos] .+= Vort3.*v1CG 

# DXa^{11}w
@views wHat .= dXdxIF[:,:,:,:,1,1].*wCG
mul!(reshape(DXwHat11,OP,OP*NF*(nz+1)),CG.DS,reshape(wHat,OP,OP*(nz+1)*NF))
# a^{21}wDY^T
@views wHat .= dXdxIF[:,:,:,:,2,1].*wCG
mul!(reshape(PermutedDimsArray(DYwHat21,(2,1,3,4)),OP,OP*NF*(nz+1)),CG.DS,reshape(PermutedDimsArray(wHat,(2,1,3,4)),OP,OP*(nz+1)*NF))
# \frac{da^{31}w}{dz}-\frac{da^{33}u}{dz}
vHat .= dXdxIC33.*v1CG
if nz>1
  @views DZuuHat31[:,:,:,1] .= 0.5 .* (vHat[:,:,:,2] .- vHat[:,:,:,1])
  @views DZuuHat31[:,:,:,2:nz] .= 0.5 .* (vHat[:,:,:,2:nz] .- vHat[:,:,:,1:nz-1])
  @views DZuuHat31[:,:,:,nz+1] .= 0.5 .* (vHat[:,:,:,nz] .- vHat[:,:,:,nz-1])
end
@views wHat .= dXdxIF[:,:,:,:,3,2] .* wCG
@views DZwwHat31 .= 0.5 .* (wHat[:,:,:,2:nz+1] .- wHat[:,:,:,1:nz])
@views FuHat[:,:,:,:,uPos] .+= 
  DZwwHat31 .* wCCG+
  0.5.*((DXwHat11[:,:,:,1:nz] .+ DYwHat21[:,:,:,1:nz] .-
  DZuuHat31[:,:,:,1:nz]).*wCG[:,:,:,1:nz] .+
  (DXwHat11[:,:,:,2:nz+1] .+ DYwHat21[:,:,:,2:nz+1] .-
  DZuuHat31[:,:,:,2:nz+1]).*wCG[:,:,:,2:nz+1])
@views FuHat[:,:,:,1:nz-1,wPos] .+=
  (-DXwHat11[:,:,:,2:nz] .- DYwHat21[:,:,:,2:nz] .+ DZuuHat31[:,:,:,2:nz]) .*
    (0.5 .* (v1CG[:,:,:,1:nz-1] .+ v1CG[:,:,:,2:nz])) .-
    0.5 .* (DZwwHat31[:,:,:,1:nz-1] .+ DZwwHat31[:,:,:,2:nz]) 

# DXa^{12}w
@views wHat .= dXdxIF[:,:,:,:,1,2].*wCG
mul!(reshape(DXwHat12,OP,OP*NF*(nz+1)),CG.DS,reshape(wHat,OP,OP*(nz+1)*NF))
# a^{22}wDY^T
@views wHat .= dXdxIF[:,:,:,:,2,2].*wCG
mul!(reshape(PermutedDimsArray(DYwHat22,(2,1,3,4)),OP,OP*NF*(nz+1)),CG.DS,reshape(PermutedDimsArray(wHat,(2,1,3,4)),OP,OP*(nz+1)*NF))
# \frac{da^{32}w}{dz}-\frac{da^{33}v}{dz}
vHat .= dXdxIC33 .* v2CG
if nz>1
  @views DZuuHat32[:,:,:,1] .= 0.5 .* (vHat[:,:,:,2] .- vHat[:,:,:,1])
  @views DZuuHat32[:,:,:,2:nz] .= 0.5 .* (vHat[:,:,:,2:nz] .- vHat[:,:,:,1:nz-1])
  @views DZuuHat32[:,:,:,nz+1] .= 0.5 .* (vHat[:,:,:,nz] .- vHat[:,:,:,nz-1])
end
@views wHat .= dXdxIF[:,:,:,:,3,2] .* wCG
@views DZwwHat32 .= 0.5 .* (wHat[:,:,:,2:nz+1] .- wHat[:,:,:,1:nz])

@views FuHat[:,:,:,:,vPos] .+= 
  (DZwwHat32).*wCCG .+
  0.5 .* ((DXwHat12[:,:,:,1:nz]  .+ DYwHat22[:,:,:,1:nz] .-
    DZuuHat32[:,:,:,1:nz]).*wCG[:,:,:,1:nz] .+
  (DXwHat12[:,:,:,2:nz+1] .+ DYwHat22[:,:,:,2:nz+1] .-
    DZuuHat32[:,:,:,2:nz+1]).*wCG[:,:,:,2:nz+1])

DZwwHat32 .= DZwwHat32 .* v2CG
@views FuHat[:,:,:,1:nz-1,wPos] .+=
    (-DXwHat12[:,:,:,2:nz] .- DYwHat22[:,:,:,2:nz] .+ DZuuHat32[:,:,:,2:nz]) .*
    (0.5 .* (v2CG[:,:,:,1:nz-1] .+ v2CG[:,:,:,2:nz])) .-
    0.5 .* (DZwwHat32[:,:,:,1:nz-1] .+ DZwwHat32[:,:,:,2:nz])
end


function FCurlNon3Vec(v1CG,v2CG,wCG,wCCG,CG,Param)
OP=CG.OrdPoly+1
NF=Param.Grid.NumFaces
nz=Param.Grid.nz
dXdxIC = Param.cache.dXdxIC
dXdxIF = Param.cache.dXdxIF
# Fu(1)
# -vv\left(DX(a^{12}u-a^{11}v)+(a^{22}u-a^{21}v)DY^T
# +\frac{d(a^{32}u-a^{31}v)}{dz}\right)
#
# ww(DXa^{11}w+a^{21}wDY^T+\frac{da^{31}w}{dz}-\frac{da^{33}u}{dz})

# Fu(2)
# uu\left(DX(a^{12}u-a^{11}v)+(a^{22}u-a^{21}v)DY^T
# +\frac{d(a^{32}u-a^{31}v)}{dz}\right)
#
# ww(DXa^{12}w+a^{22}wDY^T+\frac{da^{32}w}{dz}-\frac{d a^{33}v}{dz})

# Fw
# -uu((DXa^{11}w+a^{21}wDY^T+\frac{da^{31}w}{dz}-\frac{da^{33}u}{dz})
# -vv*((DXa^{12}w+a^{22}wDY^T+\frac{da^{32}w}{dz}-\frac{da^{33}v}{dz})

vHat1 = dXdxIC[:,:,:,:,1,2].*v1CG .- dXdxIC[:,:,:,:,1,1].*v2CG
vHat2 = dXdxIC[:,:,:,:,2,2].*v1CG .- dXdxIC[:,:,:,:,2,1].*v2CG
vHat3 = dXdxIC[:,:,:,:,3,2].*v1CG .- dXdxIC[:,:,:,:,3,1].*v2CG



DXvHat1=reshape(
  CG.DS*reshape(vHat1,OP,OP*nz*NF),
  OP,OP,NF,nz)
DYvHat2= permute(
  reshape(
  CG.DS*reshape(
  permute(vHat2
  ,[2,1,3,4])
  ,OP,OP*NF*nz)
  ,OP,OP,NF,nz)
  ,[2,1,3,4])

DZvHat3=zeros(OP,OP,NF,nz)
if nz>1
  DZvHat3[:,:,:,1]=0.5*(vHat3[:,:,:,2]-vHat3[:,:,:,1])
  DZvHat3[:,:,:,2:nz-1]=0.25*(vHat3[:,:,:,3:nz]-vHat3[:,:,:,1:nz-2])
  DZvHat3[:,:,:,nz]=0.5*(vHat3[:,:,:,nz]-vHat3[:,:,:,nz-1])
end

# DXa^{11}w
DXwHat11=reshape(CG.DS*reshape(
  dXdxIF[:,:,:,:,1,1].*wCG
  ,OP,OP*NF*(nz+1))
  ,OP,OP,NF,nz+1)
# DXa^{12}w
DXwHat12=reshape(CG.DS*reshape(
  dXdxIF[:,:,:,:,1,2].*wCG
  ,OP,OP*NF*(nz+1))
  ,OP,OP,NF,nz+1)
# a^{21}wDY^T
DYwHat21=  permute(
  reshape(
  CG.DS*reshape(
  permute(
  reshape(
  dXdxIF[:,:,:,:,2,1].*wCG
  ,OP,OP,NF,nz+1)
  ,[2,1,3,4])
  ,OP,OP*NF*(nz+1))
  ,OP,OP,NF,nz+1)
  ,[2,1,3,4])
# a^{22}wDY^T
DYwHat22= permute(
  reshape(
  CG.DS*reshape(
  permute(
  reshape(
  dXdxIF[:,:,:,:,2,2].*wCG
  ,OP,OP,NF,nz+1)
  ,[2,1,3,4])
  ,OP,OP*NF*(nz+1))
  ,OP,OP,NF,nz+1)
  ,[2,1,3,4])

# \frac{da^{31}w}{dz}-\frac{da^{33}u}{dz}
#uHat1=reshape(Param.dXdxI(:,:,:,:,3,3),OP*OP*NF,nz).*v1CG
uHat1=dXdxIC[:,:,:,:,3,3].*v1CG
DZuuHat31=zeros(OP,OP,NF,nz+1)
if nz>1
  DZuuHat31[:,:,:,1]=0.5*(uHat1[:,:,:,2]-uHat1[:,:,:,1])
  DZuuHat31[:,:,:,2:nz]=0.5*(uHat1[:,:,:,2:nz]-uHat1[:,:,:,1:nz-1])
  DZuuHat31[:,:,:,nz+1]=0.5*(uHat1[:,:,:,nz]-uHat1[:,:,:,nz-1])
end

wwHat31=dXdxIF[:,:,:,:,3,1].*wCG
DZwwHat31=0.5*(wwHat31[:,:,:,2:nz+1]-wwHat31[:,:,:,1:nz])
# \frac{da^{32}w}{dz}-\frac{da^{33}v}{dz}
#uHat2=reshape(Param.dXdxI(:,:,:,:,3,3),OP*OP*NF,nz).*v2CG
uHat2=dXdxIC[:,:,:,:,3,3].*v2CG
DZuuHat32=zeros(OP,OP,NF,nz+1)
if nz>1
  DZuuHat32[:,:,:,1]=0.5*(uHat2[:,:,:,2]-uHat2[:,:,:,1])
  DZuuHat32[:,:,:,2:nz]=0.5*(uHat2[:,:,:,2:nz]-uHat2[:,:,:,1:nz-1])
  DZuuHat32[:,:,:,nz+1]=0.5*(uHat2[:,:,:,nz]-uHat2[:,:,:,nz-1])
end
wwHat32=dXdxIF[:,:,:,:,3,2].*wCG
DZwwHat32=0.5*(wwHat32[:,:,:,2:nz+1]-wwHat32[:,:,:,1:nz])



FuHat=zeros(OP,OP,NF,nz,3)
#FuHat[:,:,:,1]=reshape(-(DXvHat1+DYvHat2+DZvHat3).*v2CG+(DXwHat11+DYwHat21+DZwHat31).*wC
#  ,OP*OP,NF,nz)
# Fu(1)
# -vv\left(DX(a^{12}u-a^{11}v)+(a^{22}u-a^{21}v)DY^T
# +\frac{d(a^{32}u-a^{31}v)}{dz}\right)
#
# ww(DXa^{11}w+a^{21}wDY^T+\frac{da^{31}w}{dz}-\frac{da^{33}u}{dz})
Vort3=DXvHat1+DYvHat2+DZvHat3
if Param.Coriolis
  str = Param.CoriolisType
  if str == "Sphere"
      lat = Param.cache.lat
      JC = Param.cache.JC
      Omega = Param.Omega
      Vort3=Vort3-reshape(2.0*Omega*sin.(
        repmat(reshape(lat[:,:,:],OP*OP*NF,1),1,nz))
        ,OP,OP,NF,nz).*
        JC
  elseif str == "Beta-Plane"
      J = Param.cache.J
      X = Param.cache.X
      beta0 = Param.beta0
      f0 = Param.f0
      y0 = Param.y0
      Vort3=Vort3-reshape((f0+beta0*(
        reshape(abs.(X[:,:,2,:,:]),OP,OP,NF,nz)-y0)).*
        J[:,:,:,:],OP*OP*NF,nz)
  end
end
FuHat[:,:,:,:,1]= -Vort3.*v2CG+
  (DZwwHat31).*wCCG+
  0.5*((DXwHat11[:,:,:,1:nz]  +DYwHat21[:,:,:,1:nz]-
  DZuuHat31[:,:,:,1:nz]).*wCG[:,:,:,1:nz]+
  (DXwHat11[:,:,:,2:nz+1]+DYwHat21[:,:,:,2:nz+1]-
  DZuuHat31[:,:,:,2:nz+1]).*wCG[:,:,:,2:nz+1])
FuHat[:,:,:,:,2]= Vort3.*v1CG+
  (DZwwHat32).*wCCG +
  0.5*((DXwHat12[:,:,:,1:nz]  +DYwHat22[:,:,:,1:nz] -
    DZuuHat32[:,:,:,1:nz]).*wCG[:,:,:,1:nz]+
  (DXwHat12[:,:,:,2:nz+1]+DYwHat22[:,:,:,2:nz+1] -
    DZuuHat32[:,:,:,2:nz+1]).*wCG[:,:,:,2:nz+1])

DZwwHat31=(DZwwHat31).*v1CG
DZwwHat32=(DZwwHat32).*v2CG
FuHat[:,:,:,1:nz-1,3]=
  (-DXwHat11[:,:,:,2:nz]-DYwHat21[:,:,:,2:nz]+DZuuHat31[:,:,:,2:nz]).*
    (0.5*(v1CG[:,:,:,1:nz-1]+v1CG[:,:,:,2:nz])) -0.5*
    (DZwwHat31[:,:,:,1:nz-1]+DZwwHat31[:,:,:,2:nz]) +
    (-DXwHat12[:,:,:,2:nz]-DYwHat22[:,:,:,2:nz]+DZuuHat32[:,:,:,2:nz]).*
    (0.5*(v2CG[:,:,:,1:nz-1]+v2CG[:,:,:,2:nz])) -0.5*
    (DZwwHat32[:,:,:,1:nz-1]+DZwwHat32[:,:,:,2:nz])

return FuHat
end
