function FDiv3LowOrderVec!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];
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

function FDiv3Vec!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];
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
  @views @. vConV = 0.5*(cCG[:,:,iz] + cCG[:,:,iz+1]) * vConV
  @views @. F[:,:,iz] -= 0.5*vConV
  @views @. F[:,:,iz+1] += 0.5*vConV
end
end

function SourceIntEnergy!(F,cCG,v1CG,v2CG,v3CG,CG,Global,iF)
@unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = Global.ThreadCache
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];
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
OP=CG.OrdPoly+1;
nz=Global.Grid.nz;
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];
@views JC = Global.Metric.JC[:,:,:,iF];
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
  @views @. vConV *= 0.5*(RhoCG[:,:,iz] + RhoCG[:,:,iz+1])
  @views @. vConV = 0.5*(abs(vConV) + vConV) * cR[:,:,iz] +
    0.5*(-abs(vConV) + vConV) * cL[:,:,iz+1];
  @views @. F[:,:,iz] -= 0.5*vConV
  @views @. F[:,:,iz+1] += 0.5*vConV
end
end

function FDiv3UpwindHorVec!(F,cCG,v1CG,v2CG,v3CG,RhoCG,CG,Global,iF)
OP=CG.OrdPoly+1;
nz=Global.Grid.nz;
@views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
@views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];
@views JC = Global.Metric.JC[:,:,:,iF];
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
    0.5*(-abs(vConV) + vConV) * cL[:,:,iz+1];
  @views @. F[:,:,iz] -= 0.5*vConV
  @views @. F[:,:,iz+1] += 0.5*vConV
end
end

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


