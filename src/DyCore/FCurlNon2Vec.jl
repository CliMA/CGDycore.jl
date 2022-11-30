function FCurlNon2Vec!(FuHat,v1CG,v2CG,CG,Global,iF)
  @unpack TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5,
  TCacheCF1, TCacheCF2, TCacheCF3, TCacheCF4 = Global.ThreadCache
  OP=CG.OrdPoly+1
  nz=Global.Grid.nz

  uPos=Global.Model.uPos
  vPos=Global.Model.vPos

  @views lat = Global.Metric.lat[:,:,iF]
  @views JC = Global.Metric.JC[:,:,:,iF] 
  @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
  @views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]

  vHat = TCacheCC1[Threads.threadid()]
  Vort3 = TCacheCC2[Threads.threadid()]
  Temp = TCacheCC3[Threads.threadid()]
#DZvHat3 = TCacheCC4[Threads.threadid()]

#wHat = TCacheCF1[Threads.threadid()]
#DXwHat11 = TCacheCF2[Threads.threadid()]
#DYwHat21 = TCacheCF3[Threads.threadid()]
#DZuuHat31 = TCacheCF4[Threads.threadid()]
#DZwwHat31 = TCacheCC5[Threads.threadid()]

#DXwHat12 = TCacheCF2[Threads.threadid()]
#DYwHat22 = TCacheCF3[Threads.threadid()]
#DZuuHat32 = TCacheCF4[Threads.threadid()]
#DZwwHat32 = TCacheCC5[Threads.threadid()]


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
end

