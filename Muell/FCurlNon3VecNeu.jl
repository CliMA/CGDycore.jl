function FCurlNon3Vec!(FuHat,v1CG,v2CG,wCG,wCCG,CG,Global,iF)
@unpack TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5,
TCacheCF1, TCacheCF2, TCacheCF3, TCacheCF4 = Global.ThreadCache
OP=CG.OrdPoly+1
nz=Global.Grid.nz

uPos=Global.Model.uPos
vPos=Global.Model.vPos
wPos = Global.Model.wPos

@views lat = Global.Metric.lat[:,:,iF]
@views JC = Global.Metric.JC[:,:,:,iF] 
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]

vHat = TCacheCC1[Threads.threadid()]
Vort3 = TCacheCC2[Threads.threadid()]
Temp = TCacheCC3[Threads.threadid()]
DZvHat3 = TCacheCC4[Threads.threadid()]

wHat = TCacheCF1[Threads.threadid()]
DXwHat11 = TCacheCF2[Threads.threadid()]
DYwHat21 = TCacheCF3[Threads.threadid()]
DZuuHat31 = TCacheCF4[Threads.threadid()]
DZwwHat31 = TCacheCC5[Threads.threadid()]

DXwHat12 = TCacheCF2[Threads.threadid()]
DYwHat22 = TCacheCF3[Threads.threadid()]
DZuuHat32 = TCacheCF4[Threads.threadid()]
DZwwHat32 = TCacheCC5[Threads.threadid()]


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

@inbounds for iz=1:nz  
  @views @. vHat[:,:,iz] = dXdxIC[:,:,iz,1,2]*v1CG[:,:,iz] - dXdxIC[:,:,iz,1,1]*v2CG[:,:,iz]
  @views mul!(Vort3[:,:,iz],CG.DS,vHat[:,:,iz])
  @views @. vHat[:,:,iz] = dXdxIC[:,:,iz,2,2]*v1CG[:,:,iz] - dXdxIC[:,:,iz,2,1]*v2CG[:,:,iz]
  @views mul!(Temp[:,:,iz],vHat[:,:,iz],CG.DST)
  @views @. Vort3[:,:,iz] += Temp[:,:,iz]
  @views @. vHat[:,:,iz] = dXdxIC[:,:,iz,3,2]*v1CG[:,:,iz] - dXdxIC[:,:,iz,3,1]*v2CG[:,:,iz]
end  
if nz>1
  @views @.Vort3[:,:,1] += 0.5*(vHat[:,:,2] - vHat[:,:,1])
  @inbounds for iz=2:nz-1
    @views @. Vort3[:,:,iz] += 0.25*(vHat[:,:,iz+1] - vHat[:,:,iz])
  end  
  @views @. Vort3[:,:,nz] += 0.5*(vHat[:,:,nz] - vHat[:,:,nz-1])
end

if Global.Model.Coriolis
  str = Global.Model.CoriolisType
  if str == "Sphere"
      Omega = Global.Phys.Omega
      @inbounds for iz = 1:nz
        @views @. Vort3[:,:,iz] -= 2.0*Omega*sin(lat)*JC[:,:,iz]  
      end   
  elseif str == "Beta-Plane"
      JC = Global.Metric.JC
      X = Global.Metric.X
      beta0 = Global.beta0
      f0 = Global.f0
      y0 = Global.y0
      @views Vort3 .= Vort3 .- (f0 .+ beta0 .* (
        reshape(abs.(0.5.*(X[:,:,1,2,:,:] .+ X[:,:,2,2,:,:])),OP,OP,NF,nz) .- y0)) .* JC
  end
end
  @views @. FuHat[:,:,:,uPos] -= Vort3[:,:,:] * v2CG[:,:,:] 
  @views @. FuHat[:,:,:,vPos] += Vort3[:,:,:] * v1CG[:,:,:] 

  @inbounds for iz=1:nz+1
    # DXa^{11}w
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,1,1] * wCG[:,:,iz]
    @views mul!(DXwHat11[:,:,iz], CG.DS, wHat[:,:,iz])
    # a^{21}wDY^T
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,2,1] * wCG[:,:,iz]
    @views mul!(DYwHat21[:,:,iz], wHat[:,:,iz], CG.DST)
  end  

  # \frac{da^{31}w}{dz}-\frac{da^{33}u}{dz}
  @views @. vHat[:,:,:] = dXdxIC[:,:,:,3,3].*v1CG[:,:,:]
  if nz>1
    @views @. DZuuHat31[:,:,1] = 0.5*(vHat[:,:,2] - vHat[:,:,1])
    @inbounds for iz=2:nz
      @views @. DZuuHat31[:,:,iz] = 0.5*(vHat[:,:,iz] - vHat[:,:,iz-1])
    end  
    @views @. DZuuHat31[:,:,nz+1] = 0.5*(vHat[:,:,nz] - vHat[:,:,nz-1])
  end
  @views @. wHat[:,:,:] = dXdxIF[:,:,:,3,1] * wCG[:,:,:]
  @inbounds for iz=1:nz
    @views @. DZwwHat31[:,:,iz] = 0.5*(wHat[:,:,iz+1] - wHat[:,:,iz])
  end  
  @views @. FuHat[:,:,:,uPos] += 
    DZwwHat31[:,:,:] * wCCG[:,:,:] +
    0.5*((DXwHat11[:,:,1:nz] + DYwHat21[:,:,1:nz] -
          DZuuHat31[:,:,1:nz]) * wCG[:,:,1:nz] +
         (DXwHat11[:,:,2:nz+1] + DYwHat21[:,:,2:nz+1] -
          DZuuHat31[:,:,2:nz+1]) * wCG[:,:,2:nz+1])

  @views @. DZwwHat31[:,:,:] = DZwwHat31[:,:,:] * v1CG[:,:,:]
  @views @. FuHat[:,:,1:nz-1,wPos] +=
    (-DXwHat11[:,:,2:nz] - DYwHat21[:,:,2:nz] + DZuuHat31[:,:,2:nz]) *
    (0.5*(v1CG[:,:,1:nz-1] + v1CG[:,:,2:nz])) -
    0.5*(DZwwHat31[:,:,1:nz-1]+DZwwHat31[:,:,2:nz])

  @inbounds for iz=1:nz+1
    # DXa^{12}w
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,1,2] * wCG[:,:,iz]
    @views mul!(DXwHat12[:,:,iz], CG.DS, wHat[:,:,iz])
    # a^{22}wDY^T
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,2,2] * wCG[:,:,iz]
    @views mul!(DYwHat22[:,:,iz], wHat[:,:,iz], CG.DST)
  end  

  # \frac{da^{32}w}{dz}-\frac{da^{33}v}{dz}
  @views @. vHat[:,:,:] = dXdxIC[:,:,:,3,3] * v2CG[:,:,:]
  if nz>1
    @views @. DZuuHat32[:,:,1] = 0.5*(vHat[:,:,2] - vHat[:,:,1])
    @inbounds for iz=2:nz
      @views @. DZuuHat32[:,:,iz] = 0.5*(vHat[:,:,iz] - vHat[:,:,iz-1])
    end
    @views @. DZuuHat32[:,:,nz+1] = 0.5*(vHat[:,:,nz] - vHat[:,:,nz-1])
  end
  @views @. wHat[:,:,:] = dXdxIF[:,:,:,3,2] * wCG[:,:,:]
  @views @. DZwwHat32[:,:,:] = 0.5*(wHat[:,:,2:nz+1] - wHat[:,:,1:nz])

  @views @. FuHat[:,:,:,vPos] += 
    DZwwHat32[:,:,:] * wCCG[:,:,:] +
    0.5*((DXwHat12[:,:,1:nz]  + DYwHat22[:,:,1:nz] -
          DZuuHat32[:,:,1:nz]) * wCG[:,:,1:nz] +
         (DXwHat12[:,:,2:nz+1] + DYwHat22[:,:,2:nz+1] -
          DZuuHat32[:,:,2:nz+1]) * wCG[:,:,2:nz+1])

  @views @. DZwwHat32[:,:,:] = DZwwHat32[:,:,:] * v2CG[:,:,:]
  @views @. FuHat[:,:,1:nz-1,wPos] +=
    (-DXwHat12[:,:,2:nz] - DYwHat22[:,:,2:nz] + DZuuHat32[:,:,2:nz]) *
    (0.5*(v2CG[:,:,1:nz-1] + v2CG[:,:,2:nz])) -
    0.5*(DZwwHat32[:,:,1:nz-1] + DZwwHat32[:,:,2:nz])
end

function FCurlNon3VecTest!(FuHat1,v1CG,v2CG,wCG,wCCG,CG,Global,iF)
@unpack TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5,
TCacheCF1, TCacheCF2, TCacheCF3, TCacheCF4 = Global.ThreadCache
OP=CG.OrdPoly+1
nz=Global.Grid.nz
FuHat = similar(FuHat1)
FuHat .= 0.0

uPos=Global.Model.uPos
vPos=Global.Model.vPos
wPos = Global.Model.wPos

@views lat = Global.Metric.lat[:,:,iF]
@views JC = Global.Metric.JC[:,:,:,iF] 
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]

vHat = TCacheCC1[Threads.threadid()]
Vort3 = TCacheCC2[Threads.threadid()]
Temp = TCacheCC3[Threads.threadid()]
DZvHat3 = TCacheCC4[Threads.threadid()]

wHat = TCacheCF1[Threads.threadid()]
DXwHat11 = TCacheCF2[Threads.threadid()]
DYwHat21 = TCacheCF3[Threads.threadid()]
DZuuHat31 = TCacheCF4[Threads.threadid()]
DZwwHat31 = TCacheCC5[Threads.threadid()]

DXwHat12 = TCacheCF2[Threads.threadid()]
DYwHat22 = TCacheCF3[Threads.threadid()]
DZuuHat32 = TCacheCF4[Threads.threadid()]
DZwwHat32 = TCacheCC5[Threads.threadid()]


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

@inbounds for iz=1:nz  
  @views @. vHat[:,:,iz] = dXdxIC[:,:,iz,1,2]*v1CG[:,:,iz] - dXdxIC[:,:,iz,1,1]*v2CG[:,:,iz]
  @views mul!(Vort3[:,:,iz],CG.DS,vHat[:,:,iz])
  @views @. vHat[:,:,iz] = dXdxIC[:,:,iz,2,2]*v1CG[:,:,iz] - dXdxIC[:,:,iz,2,1]*v2CG[:,:,iz]
  @views mul!(Temp[:,:,iz],vHat[:,:,iz],CG.DST)
  @views @. Vort3[:,:,iz] += Temp[:,:,iz]
  @views @. vHat[:,:,iz] = dXdxIC[:,:,iz,3,2]*v1CG[:,:,iz] - dXdxIC[:,:,iz,3,1]*v2CG[:,:,iz]
end  
if nz>1
  @views @. Vort3[:,:,1] += 0.5*(vHat[:,:,2] - vHat[:,:,1])
  @inbounds for iz=2:nz-1
    @views @. Vort3[:,:,iz] += 0.25*(vHat[:,:,iz+1] - vHat[:,:,iz])
  end  
  @views @. Vort3[:,:,nz] += 0.5*(vHat[:,:,nz] - vHat[:,:,nz-1])
end

if Global.Model.Coriolis
  str = Global.Model.CoriolisType
  if str == "Sphere"
      Omega = Global.Phys.Omega
      @inbounds for iz = 1:nz
        @views @. Vort3[:,:,iz] -= 2.0*Omega*sin(lat)*JC[:,:,iz]  
      end   
  elseif str == "Beta-Plane"
      JC = Global.Metric.JC
      X = Global.Metric.X
      beta0 = Global.beta0
      f0 = Global.f0
      y0 = Global.y0
      @views Vort3 .= Vort3 .- (f0 .+ beta0 .* (
        reshape(abs.(0.5.*(X[:,:,1,2,:,:] .+ X[:,:,2,2,:,:])),OP,OP,NF,nz) .- y0)) .* JC
  end
end
  @views @. FuHat[:,:,:,uPos] -= Vort3[:,:,:] * v2CG[:,:,:] 
  @views @. FuHat[:,:,:,vPos] += Vort3[:,:,:] * v1CG[:,:,:] 

  @inbounds for iz=1:nz+1
    # DXa^{11}w
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,1,1] * wCG[:,:,iz]
    @views mul!(DXwHat11[:,:,iz], CG.DS, wHat[:,:,iz])
    # a^{21}wDY^T
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,2,1] * wCG[:,:,iz]
    @views mul!(DYwHat21[:,:,iz], wHat[:,:,iz], CG.DST)
  end  

  # \frac{da^{31}w}{dz}-\frac{da^{33}u}{dz}
  @views @. vHat[:,:,:] = dXdxIC[:,:,:,3,3].*v1CG[:,:,:]
  if nz>1
    @views @. DZuuHat31[:,:,1] = 0.5*(vHat[:,:,2] - vHat[:,:,1])
    @inbounds for iz=2:nz
      @views @. DZuuHat31[:,:,iz] = 0.5*(vHat[:,:,iz] - vHat[:,:,iz-1])
    end  
    @views @. DZuuHat31[:,:,nz+1] = 0.5*(vHat[:,:,nz] - vHat[:,:,nz-1])
  end
  @views @. wHat[:,:,:] = dXdxIF[:,:,:,3,1] * wCG[:,:,:]
  @inbounds for iz=1:nz
    @views @. DZwwHat31[:,:,iz] = 0.5*(wHat[:,:,iz+1] - wHat[:,:,iz])
  end  
  @views @. FuHat[:,:,:,uPos] += 
    DZwwHat31[:,:,:] * wCCG[:,:,:] +
    0.5*((DXwHat11[:,:,1:nz] + DYwHat21[:,:,1:nz] -
          DZuuHat31[:,:,1:nz]) * wCG[:,:,1:nz] +
         (DXwHat11[:,:,2:nz+1] + DYwHat21[:,:,2:nz+1] -
          DZuuHat31[:,:,2:nz+1]) * wCG[:,:,2:nz+1])

  @views @. DZwwHat31[:,:,:] = DZwwHat31[:,:,:] * v1CG[:,:,:]
  @views @. FuHat[:,:,1:nz-1,wPos] +=
    (-DXwHat11[:,:,2:nz] - DYwHat21[:,:,2:nz] + DZuuHat31[:,:,2:nz]) *
    (0.5*(v1CG[:,:,1:nz-1] + v1CG[:,:,2:nz])) -
    0.5*(DZwwHat31[:,:,1:nz-1]+DZwwHat31[:,:,2:nz])

  @inbounds for iz=1:nz+1
    # DXa^{12}w
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,1,2] * wCG[:,:,iz]
    @views mul!(DXwHat12[:,:,iz], CG.DS, wHat[:,:,iz])
    # a^{22}wDY^T
    @views @. wHat[:,:,iz] = dXdxIF[:,:,iz,2,2] * wCG[:,:,iz]
    @views mul!(DYwHat22[:,:,iz], wHat[:,:,iz], CG.DST)
  end  

  # \frac{da^{32}w}{dz}-\frac{da^{33}v}{dz}
  @views @. vHat[:,:,:] = dXdxIC[:,:,:,3,3] * v2CG[:,:,:]
  if nz>1
    @views @. DZuuHat32[:,:,1] = 0.5*(vHat[:,:,2] - vHat[:,:,1])
    @inbounds for iz=2:nz
      @views @. DZuuHat32[:,:,iz] = 0.5*(vHat[:,:,iz] - vHat[:,:,iz-1])
    end
    @views @. DZuuHat32[:,:,nz+1] = 0.5*(vHat[:,:,nz] - vHat[:,:,nz-1])
  end
  @views @. wHat[:,:,:] = dXdxIF[:,:,:,3,2] * wCG[:,:,:]
  @views @. DZwwHat32[:,:,:] = 0.5*(wHat[:,:,2:nz+1] - wHat[:,:,1:nz])

  @views @. FuHat[:,:,:,vPos] += 
    DZwwHat32[:,:,:] * wCCG[:,:,:] +
    0.5*((DXwHat12[:,:,1:nz]  + DYwHat22[:,:,1:nz] -
          DZuuHat32[:,:,1:nz]) * wCG[:,:,1:nz] +
         (DXwHat12[:,:,2:nz+1] + DYwHat22[:,:,2:nz+1] -
          DZuuHat32[:,:,2:nz+1]) * wCG[:,:,2:nz+1])

  @views @. DZwwHat32[:,:,:] = DZwwHat32[:,:,:] * v2CG[:,:,:]
  @views @. FuHat[:,:,1:nz-1,wPos] +=
    (-DXwHat12[:,:,2:nz] - DYwHat22[:,:,2:nz] + DZuuHat32[:,:,2:nz]) *
    (0.5*(v2CG[:,:,1:nz-1] + v2CG[:,:,2:nz])) -
    0.5*(DZwwHat32[:,:,1:nz-1] + DZwwHat32[:,:,2:nz])
  s=sum(v1CG .* FuHat[:,:,:,uPos] .+ v2CG .* FuHat[:,:,:,vPos]) +sum( wCG[:,:,2:nz] .* FuHat[:,:,1:nz-1,wPos])   
  @show s
  s1=sum(abs.(v1CG .* FuHat[:,:,:,uPos]) .+ abs.(v2CG .* FuHat[:,:,:,vPos])) +
     sum( abs.(wCG[:,:,2:nz] .* FuHat[:,:,1:nz-1,wPos]))   
  @show s1
end

