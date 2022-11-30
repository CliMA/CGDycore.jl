function FDiv3LowOrderVec!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1
NF=Global.Grid.NumFaces
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
# Contravariant components

vCon1 = TCacheC1[Threads.threadid()]
vCon2 = TCacheC2[Threads.threadid()]
DvCon1 = TCacheC3[Threads.threadid()]
DvCon2 = TCacheC4[Threads.threadid()]
vConV = TCacheC1[Threads.threadid()]

@inbounds for iz=1:nz
  @views @. vCon1 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) 
  @views @. vCon1[1:OP-1,:] = 0.5 * (vCon1[1:OP-1,:] + vCon1[2:OP,:])
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) 
  @views @. vCon2[:,OP-1] = 0.5 * (vCon2[:,1:OP-1] + vCon2[:,2:OP])
  @. DvCon1 = 0.0
  @views @. DvCon1[1:OP-1,:] += -0.5 * (abs(vCon1[1:OP-1,:])+vCon1[1:OP-1,:])*cCG[1:OP-1,:] 
  @views @. DvCon1[2:OP,:] += 0.5 * (abs(vCon1[1:OP-1,:])-vCon1[1:OP-1,:])*cCG[2:OP,:] 
  @. DvCon2 = 0.0
  @views @. DvCon2[:,1:OP-1] += -0.5 * (abs(vCon1[:,1:OP-1])+vCon1[:,1:OP-1])*cCG[:,1:OP-1] 
  @views @. DvCon2[:,2:OP] += 0.5 * (abs(vCon2[:,1:OP-1])-vCon2[:,1:OP-1])*cCG[:,2:OP] 
  @views @. F[:,:,iz]  = F[:,:,iz]  - DvCon1 - DvCon2 
end
@inbounds for iz=1:nz-1
  @views @. vConV = 0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
  @views @. vConV = 0.5*(cCG[:,:,iz] + cCG[:,:,iz+1]) * vConV
  @views @. F[:,:,iz] -= 0.5*vConV
  @views @. F[:,:,iz+1] += 0.5*vConV
end
end

function AdvecVertVec!(f,RhoCG,v1CG,v2CG,v3CG,Omega,CG,Global,iF)
  nz=Global.Grid.nz
  @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
  @views JC = Global.Metric.JC[:,:,:,iF]
  @inbounds for iz=1:nz-1
    @views @. Omega[:,:,iz+1] = (v1CG[:,:,iz]*dXdxIC[:,:,iz,3,1] + v1CG[:,:,iz+1]*dXdxIC[:,:,iz+1,3,1] +
                               v2CG[:,:,iz]*dXdxIC[:,:,iz,3,2] + v2CG[:,:,iz+1]*dXdxIC[:,:,iz+1,3,2] +
       v3CG[:,:,iz+1] ) /  (dXdxIC[:,:,iz,3,3] + dXdxIC[:,:,iz,3,3]) * 
       (RhoCG[:,:,iz] * JC[:,:,iz] + RhoCG[:,:,iz+1] * JC[:,:,iz+1]) / (JC[:,:,iz] + JC[:,:,iz+1])
  end   
  @views @. dz = JC[:,:,:,iF] / dXdxIC[:,:,:,3,3,iF]
  iz = 1
  @views @. f[:,:,iz,1] += JC[:,:,iz]/RhoCG[:,:,iz] *
    ((Omega[:,:,iz+1]*(v1CG[:,:,iz+1]+v1CG[:,:,iz]) -
    Omega[:,:,iz]*(v1CG[:,:,iz]+v1CG[:,:,iz])) /(2.0*dz[:,:,iz]) -
    v1CG[:,:,iz]*(Omega[:,:,iz+1]-Omega[:,:,iz]) / dz[:,:,iz]);
  @views @. f[:,:,iz,2] += JC[:,:,iz]/RhoCG[:,:,iz] *
    ((Omega[:,:,iz+1]*(v2CG[:,:,iz+1]+v2CG[:,:,iz]) -
    Omega[:,:,iz]*(v2CG[:,:,iz]+v2CG[:,:,iz])) /(2.0*dz[:,:,iz]) -
    v1CG[:,:,iz]*(Omega[:,:,iz+1]-Omega[:,:,iz]) / dz[:,:,iz]);
  for iz = 2 : nz-1
    @views @. f[:,:,iz,1] += JC[:,:,iz]/RhoCG[:,:,iz] *
      ((Omega[:,:,iz+1]*(v1CG[:,:,iz+1]+v1CG[:,:,iz]) -
      Omega[:,:,iz]*(v1CG[:,:,iz]+v1CG[:,:,iz-1])) /(2.0*dz[:,:,iz]) -
      v1CG[:,:,iz]*(Omega[:,:,iz+1]-Omega[:,:,iz]) / dz[:,:,iz]);
    @views @. f[:,:,iz,2] += JC[:,:,iz]/RhoCG[:,:,iz] *
      ((Omega[:,:,iz+1]*(v2CG[:,:,iz+1]+v2CG[:,:,iz]) -
      Omega[:,:,iz]*(v2CG[:,:,iz]+v2CG[:,:,iz-1])) /(2.0*dz[:,:,iz]) -
      v2CG[:,:,iz]*(Omega[:,:,iz+1]-Omega[:,:,iz]) / dz[:,:,iz]);
  end
  iz = nz
  @views @. f[:,:,iz,1] += JC[:,:,iz]/RhoCG[:,:,iz] *
    ((Omega[:,:,iz+1]*(v1CG[:,:,iz]+v1CG[:,:,iz]) -
    Omega[:,:,iz]*(v1CG[:,:,iz]+v1CG[:,:,iz-1])) /(2.0*dz[:,:,iz]) -
    v1CG[:,:,iz]*(Omega[:,:,iz+1]-Omega[:,:,iz]) / dz[:,:,iz]);
  @views @. f[:,:,iz,2] += JC[:,:,iz]/RhoCG[:,:,iz] *
    ((Omega[:,:,iz+1]*(v2CG[:,:,iz]+v2CG[:,:,iz]) -
    Omega[:,:,iz]*(v2CG[:,:,iz]+v2CG[:,:,iz-1])) /(2.0*dz[:,:,iz]) -
    v2CG[:,:,iz]*(Omega[:,:,iz+1]-Omega[:,:,iz]) / dz[:,:,iz]);
end  


function FDiv3Vec!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1
NF=Global.Grid.NumFaces
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
@views JC = Global.Metric.JC[:,:,:,iF]
# Contravariant components

vCon1 = TCacheC1[Threads.threadid()]
vCon2 = TCacheC2[Threads.threadid()]
DvCon1 = TCacheC3[Threads.threadid()]
DvCon2 = TCacheC4[Threads.threadid()]
vConV = TCacheC1[Threads.threadid()]

@inbounds for iz=1:nz
  @views @. vCon1 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) * cCG[:,:,iz]
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) * cCG[:,:,iz]
  mul!(DvCon1,CG.DS,vCon1)
  mul!(DvCon2,vCon2,CG.DST)
  @views @. F[:,:,iz]  = F[:,:,iz]  - DvCon1 - DvCon2 
end
@inbounds for iz=1:nz-1
  @views @. vConV = 0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
# @views @. vConV = 0.5*(cCG[:,:,iz] + cCG[:,:,iz+1]) * vConV
  @views @. vConV *= (cCG[:,:,iz] * JC[:,:,iz] + cCG[:,:,iz+1] * JC[:,:,iz+1]) / (JC[:,:,iz] + JC[:,:,iz+1]) 
  @views @. F[:,:,iz] = F[:,:,iz] - 0.5*vConV
  @views @. F[:,:,iz+1] = F[:,:,iz+1] + 0.5*vConV
end
end

function FDiv3ExpVec!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1
NF=Global.Grid.NumFaces
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
# Contravariant components

vCon1 = TCacheC1[Threads.threadid()]
vCon2 = TCacheC2[Threads.threadid()]
DvCon1 = TCacheC3[Threads.threadid()]
DvCon2 = TCacheC4[Threads.threadid()]
vConV = TCacheC1[Threads.threadid()]

@inbounds for iz=1:nz
  @views @. vCon1 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) * cCG[:,:,iz]
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) * cCG[:,:,iz]
  mul!(DvCon1,CG.DS,vCon1)
  mul!(DvCon2,vCon2,CG.DST)
  @views @. F[:,:,iz]  = F[:,:,iz]  - DvCon1 - DvCon2 
end
@inbounds for iz=1:nz-1
  @views @. vConV = 0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) 
#    v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
  @views @. vConV = 0.5*(cCG[:,:,iz] + cCG[:,:,iz+1]) * vConV
  @views @. F[:,:,iz] = F[:,:,iz] - 0.5*vConV
  @views @. F[:,:,iz+1] = F[:,:,iz+1] + 0.5*vConV
end
end

function FDiv3ImpVec!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1
NF=Global.Grid.NumFaces
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
# Contravariant components

vCon1 = TCacheC1[Threads.threadid()]
vCon2 = TCacheC2[Threads.threadid()]
DvCon1 = TCacheC3[Threads.threadid()]
DvCon2 = TCacheC4[Threads.threadid()]
vConV = TCacheC1[Threads.threadid()]

@inbounds for iz=1:nz-1
  @views @. vConV = #0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
#    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
  @views @. vConV = 0.5*(cCG[:,:,iz] + cCG[:,:,iz+1]) * vConV
  @views @. F[:,:,iz] = F[:,:,iz] - 0.5*vConV
  @views @. F[:,:,iz+1] = F[:,:,iz+1] + 0.5*vConV
end
end

function FDiv3ImpGlobalVec!(F,c,v3,Global,iG)
  nz=Global.Grid.nz
  @views dz = Global.Metric.dz[:,iG]

  @inbounds for iz=1:nz-1
    vConV = 0.5*(c[iz] + c[iz+1]) * v3[iz]
    F[iz] = F[iz] - vConV / dz[iz]
    F[iz+1] = F[iz+1] + vConV / dz[iz+1]
  end
end

function SourceIntEnergy!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1
NF=Global.Grid.NumFaces
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
# Contravariant components

vCon1 = TCacheC1[Threads.threadid()]
vCon2 = TCacheC2[Threads.threadid()]
DvCon1 = TCacheC3[Threads.threadid()]
DvCon2 = TCacheC4[Threads.threadid()]
vConV = TCacheC1[Threads.threadid()]

@inbounds for iz=1:nz
  @views @. vCon1 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) 
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) 
  mul!(DvCon1,CG.DS,vCon1)
  mul!(DvCon2,vCon2,CG.DST)
  @views @. F[:,:,iz]  = F[:,:,iz] - (DvCon1 + DvCon2) * cCG[:,:,iz] 
end
@inbounds for iz=1:nz-1
  @views @. vConV = 0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
  @views @. F[:,:,iz] = F[:,:,iz] - 0.5*vConV * cCG[:,:,iz]
  @views @. F[:,:,iz+1] = F[:,:,iz+1] + 0.5*vConV * cCG[:,:,iz+1]
end
end

function FDiv3UpwindVec!(F,cCG,v1CG,v2CG,v3CG,RhoCG,CG,Global,iF)
OP=CG.OrdPoly+1
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
@views JC = Global.Metric.JC[:,:,:,iF]
# Contravariant components

vCon1 = Global.Cache.CacheE1
vCon2 = Global.Cache.CacheE2
DvCon1 = Global.Cache.CacheE3
DvCon2 = Global.Cache.CacheE4
vConV = Global.Cache.CacheE1
cL = Global.Cache.CacheC1
cR = Global.Cache.CacheC2
qCG = Global.Cache.CacheC3

@inbounds for iz=1:nz
  @views @. vCon1 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) * cCG[:,:,iz]
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) * cCG[:,:,iz]
  mul!(DvCon1,CG.DS,vCon1)
  mul!(DvCon2,vCon2,CG.DST)
  @views @. F[:,:,iz]  = F[:,:,iz]  - DvCon1 - DvCon2 
end

@. qCG = cCG / RhoCG
if nz>1
  @inbounds for j=1:OP
    @inbounds for i=1:OP  
      qCG0 = ((3.0 * qCG[i,j,1] - 2.0 * qCG[i,j,2]) * JC[i,j,1] + qCG[i,j,1] * JC[i,j,2]) / (JC[i,j,1] + JC[i,j,2]) 
      (cL[i,j,1],cR[i,j,1]) = Rec3(qCG0,qCG[i,j,1],qCG[i,j,2],JC[i,j,1],JC[i,j,1],JC[i,j,2])  
    end
  end  
  @inbounds for iz=2:nz-1
    @inbounds for j=1:OP
      @inbounds for i=1:OP  
        (cL[i,j,iz],cR[i,j,iz]) = Rec3(qCG[i,j,iz-1],qCG[i,j,iz],qCG[i,j,iz+1],JC[i,j,iz-1],JC[i,j,iz],JC[i,j,iz+1])  
      end
    end  
  end  
  @inbounds for j=1:OP
    @inbounds for i=1:OP  
      qCG1 = ((3.0 * qCG[i,j,nz] - 2.0 * qCG[i,j,nz-1]) * JC[i,j,nz] + qCG[i,j,nz] * JC[i,j,nz-1]) / (JC[i,j,nz-1] + JC[i,j,nz]) 
      (cL[i,j,nz],cR[i,j,nz]) = Rec3(qCG[i,j,nz-1],qCG[i,j,nz],qCG1,JC[i,j,nz-1],JC[i,j,nz],JC[i,j,nz])  
    end
  end  
end
@inbounds for iz=1:nz-1
  @views @. vConV = 0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
# @views @. vConV *= 0.5*(RhoCG[:,:,iz] + RhoCG[:,:,iz+1])
  @views @. vConV *= (RhoCG[:,:,iz] * JC[:,:,iz] + RhoCG[:,:,iz+1] * JC[:,:,iz+1]) / (JC[:,:,iz] + JC[:,:,iz+1]) 
  @views @. vConV = 0.5*(abs(vConV) + vConV) * cR[:,:,iz] +
    0.5*(-abs(vConV) + vConV) * cL[:,:,iz+1]
  @views @. F[:,:,iz] = F[:,:,iz] - 0.5*vConV
  @views @. F[:,:,iz+1] = F[:,:,iz+1] + 0.5*vConV
end
end

function FDiv3UpwindImpVec!(F,cCG,v1CG,v2CG,v3CG,RhoCG,CG,Global,iF)
OP=CG.OrdPoly+1
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
@views JC = Global.Metric.JC[:,:,:,iF]
# Contravariant components

vCon1 = Global.Cache.CacheE1
vCon2 = Global.Cache.CacheE2
DvCon1 = Global.Cache.CacheE3
DvCon2 = Global.Cache.CacheE4
vConV = Global.Cache.CacheE1
cL = Global.Cache.CacheC1
cR = Global.Cache.CacheC2
qCG = Global.Cache.CacheC3


@. qCG = cCG / RhoCG
if nz>1
  @inbounds for j=1:OP
    @inbounds for i=1:OP  
      qCG0 = ((3.0 * qCG[i,j,1] - 2.0 * qCG[i,j,2]) * JC[i,j,1] + qCG[i,j,1] * JC[i,j,2]) / (JC[i,j,1] + JC[i,j,2]) 
      (cL[i,j,1],cR[i,j,1]) = Rec3(qCG0,qCG[i,j,1],qCG[i,j,2],JC[i,j,1],JC[i,j,1],JC[i,j,2])  
    end
  end  
  @inbounds for iz=2:nz-1
    @inbounds for j=1:OP
      @inbounds for i=1:OP  
        (cL[i,j,iz],cR[i,j,iz]) = Rec3(qCG[i,j,iz-1],qCG[i,j,iz],qCG[i,j,iz+1],JC[i,j,iz-1],JC[i,j,iz],JC[i,j,iz+1])  
      end
    end  
  end  
  @inbounds for j=1:OP
    @inbounds for i=1:OP  
      qCG1 = ((3.0 * qCG[i,j,nz] - 2.0 * qCG[i,j,nz-1]) * JC[i,j,nz] + qCG[i,j,nz] * JC[i,j,nz-1]) / (JC[i,j,nz-1] + JC[i,j,nz]) 
      (cL[i,j,nz],cR[i,j,nz]) = Rec3(qCG[i,j,nz-1],qCG[i,j,nz],qCG1,JC[i,j,nz-1],JC[i,j,nz],JC[i,j,nz])  
    end
  end  
end
@inbounds for iz=1:nz-1
  @views @. vConV = #0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
#    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
# @views @. vConV *= 0.5*(RhoCG[:,:,iz] + RhoCG[:,:,iz+1])
  @views @. vConV *= (RhoCG[:,:,iz] * JC[:,:,iz] + RhoCG[:,:,iz+1] * JC[:,:,iz+1]) / (JC[:,:,iz] + JC[:,:,iz+1]) 
  @views @. vConV = 0.5*(abs(vConV) + vConV) * cR[:,:,iz] +
    0.5*(-abs(vConV) + vConV) * cL[:,:,iz+1]
  @views @. F[:,:,iz] = F[:,:,iz] - 0.5*vConV
  @views @. F[:,:,iz+1] = F[:,:,iz+1] + 0.5*vConV
end
end

function FDiv3UpwindImpGlobalVec!(F,c,v3,Rho,Global,iG)
  nz=Global.Grid.nz
  @views dz = Global.Metric.dz[:,iG]
  @views cL = Global.Cache.CacheC1[1,1,:]
  @views cR = Global.Cache.CacheC2[1,1,:]
  @views q = Global.Cache.CacheC3[1,1,:]


  @. q = c / Rho
  if nz>1
    q0 = ((3.0 * q[1] - 2.0 * q[2]) * dz[1] + q[1] * dz[2]) / (dz[1] + dz[2]) 
    (cL[1],cR[1]) = Rec3(q0,q[1],q[2],dz[1],dz[1],dz[2])  
    @inbounds for iz=2:nz-1
      (cL[iz],cR[iz]) = Rec3(q[iz-1],q[iz],q[iz+1],dz[iz-1],dz[iz],dz[iz+1])  
    end  
    q1 = ((3.0 * q[nz] - 2.0 * q[nz-1]) * dz[nz] + q[nz] * dz[nz-1]) / (dz[nz-1] + dz[nz]) 
    (cL[nz],cR[nz]) = Rec3(q[nz-1],q[nz],q1,dz[nz-1],dz[nz],dz[nz])  
  end  
  @inbounds for iz=1:nz-1
    vConV = (Rho[iz] * dz[iz] + Rho[iz+1] * dz[iz+1]) / (dz[iz] + dz[iz1]) * v3[iz]
    vConV = 0.5*(abs(vConV) + vConV) * cR[iz] +
    0.5*(-abs(vConV) + vConV) * cL[iz+1]
    F[iz] = F[iz] - vConV / dz[iz]
    F[iz+1] = F[iz+1] + vConV / dz[iz+1]
  end
end

function FDiv3UpwindHorVec!(F,cCG,v1CG,v2CG,v3CG,RhoCG,CG,Global,iF)
OP=CG.OrdPoly+1
nz=Global.Grid.nz
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF]
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF]
@views JC = Global.Metric.JC[:,:,:,iF]
# Contravariant components

vCon = Global.Cache.CacheE1
DvCon = Global.Cache.CacheE2
DvC = Global.Cache.CacheE3
DvRho = Global.Cache.CacheE4
vConV = Global.Cache.CacheE1
cL = Global.Cache.CacheC1
cR = Global.Cache.CacheC2
qCG = Global.Cache.CacheC3

@inbounds for iz=1:nz
  @views @. vCon = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) * cCG[:,:,iz]
  mul!(DvC,CG.DS,vCon1)
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) * cCG[:,:,iz]
  mul!(DvCon,vCon2,CG.DST)
  @views @. Dvc += DvCon
  @views @. vCon = (v1CG[:,:,iz] * dXdxIC[:,:,iz,1,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,1,2]) 
  mul!(DvRho,CG.DS,vCon1)
  @views @. vCon2 = (v1CG[:,:,iz] * dXdxIC[:,:,iz,2,1] + 
    v2CG[:,:,iz] * dXdxIC[:,:,iz,2,2]) 
  mul!(DvCon,vCon2,CG.DST)
  @views @. DvRho += DvCon
  
end
  #@views @. F[:,:,iz]  = F[:,:,iz]  - DvCon1 - DvCon2 

@. qCG = cCG / RhoCG
if nz>1
  @inbounds for j=1:OP
    @inbounds for i=1:OP  
      (cL[i,j,1],cR[i,j,1]) = Rec3(qCG[i,j,1],qCG[i,j,1],qCG[i,j,2],JC[i,j,1],JC[i,j,1],JC[i,j,2])  
    end
  end  
  @inbounds for iz=2:nz-1
    @inbounds for j=1:OP
      @inbounds for i=1:OP  
        (cL[i,j,iz],cR[i,j,iz]) = Rec3(qCG[i,j,iz-1],qCG[i,j,iz],qCG[i,j,iz+1],JC[i,j,iz-1],JC[i,j,iz],JC[i,j,iz+1])  
      end
    end  
  end  
  @inbounds for j=1:OP
    @inbounds for i=1:OP  
      (cL[i,j,nz],cR[i,j,nz]) = Rec3(qCG[i,j,nz-1],qCG[i,j,nz],qCG[i,j,nz],JC[i,j,nz-1],JC[i,j,nz],JC[i,j,nz])  
    end
  end  
end
@inbounds for iz=1:nz-1
  @views @. vConV = 0.5*((v1CG[:,:,iz] + v1CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,1] +
    (v2CG[:,:,iz] + v2CG[:,:,iz+1]) * dXdxIF[:,:,iz+1,3,2]) +
     v3CG[:,:,iz+1] * dXdxIF[:,:,iz+1,3,3]
  @views @. vConV *= 0.5*(RhoCG[:,:,iz] + RhoCG[:,:,iz+1])
  @views @. vConV = 0.5*(abs(vConV) + vConV) * cR[:,:,iz] +
    0.5*(-abs(vConV) + vConV) * cL[:,:,iz+1]
  @views @. F[:,:,iz] -= 0.5*vConV
  @views @. F[:,:,iz+1] += 0.5*vConV
end
end

#function Rec3(cL,cC,cR,JL,JC,JR)
function Rec3(cL,cC,cR,JL,JC,JR)
  kL=(JC/(JC+JL))*((JR+JC)/(JL+JC+JR))
  kR=-(JC/(JR+JC))*(JL/(JL+JC+JR))
  cCL=kL*cL+(1.0-kL-kR)*cC+kR*cR
  kR=(JC/(JC+JR))*((JL+JC)/(JL+JC+JR))
  kL=-(JC/(JL+JC))*(JR/(JL+JC+JR))
  cCR=kL*cL+(1.0-kL-kR)*cC+kR*cR
  return (cCL,cCR)
# return(cC,cC)
end
function Rec1(cL,cC,cR,JL,JC,JR)
  return(cC,cC)
end


