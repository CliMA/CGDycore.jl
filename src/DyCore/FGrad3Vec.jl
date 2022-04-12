function FGrad3Vec!(F,cCG,CG,Global,iF)
    (;  uPos,
        vPos,
        wPos) = Global.Model
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;

  D1cCG = Global.Cache.CacheC1
  D2cCG = Global.Cache.CacheC2
  D3cCG = Global.Cache.CacheC3
  D3cCGE = Global.Cache.CacheC4
  @views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
  @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];


  @inbounds for iz=1:nz-1
    @views @. D3cCG[:,:,iz] = 0.5*(cCG[:,:,iz+1] - cCG[:,:,iz])
    @views @. F[:,:,iz,wPos] -= dXdxIF[:,:,iz+1,3,3]*D3cCG[:,:,iz] 
  end  
  @views mul!(D1cCG[:,:,1],CG.DS,cCG[:,:,1])
  @views mul!(D2cCG[:,:,1],cCG[:,:,1],CG.DST)
  if nz>1
    @views  @. D3cCGE[:,:,1] = 1.5*D3cCG[:,:,1] - 0.5*D3cCG[:,:,2]
    @views mul!(D1cCG[:,:,nz],CG.DS,cCG[:,:,nz])
    @views mul!(D2cCG[:,:,nz],cCG[:,:,nz],CG.DST)
    @views  @. D3cCGE[:,:,nz] = D3cCG[:,:,nz-1]
    @views @. F[:,:,nz,uPos] -= (dXdxIC[:,:,nz,1,1]*D1cCG[:,:,nz] +
      dXdxIC[:,:,nz,2,1]*D2cCG[:,:,nz] +
      dXdxIC[:,:,nz,3,1]*D3cCGE[:,:,nz] )
    @views @. F[:,:,nz,vPos] -= (dXdxIC[:,:,nz,1,2]*D1cCG[:,:,nz] +
      dXdxIC[:,:,nz,2,2]*D2cCG[:,:,nz] +
      dXdxIC[:,:,nz,3,2]*D3cCGE[:,:,nz] )
  else
    @views  @. D3cCGE[:,:,1] = 0.0  
  end  
  @views @. F[:,:,1,uPos] -= (dXdxIC[:,:,1,1,1]*D1cCG[:,:,1] +
    dXdxIC[:,:,1,2,1]*D2cCG[:,:,1] +
    dXdxIC[:,:,1,3,1]*D3cCGE[:,:,1] ) 
  @views @. F[:,:,1,vPos] -= (dXdxIC[:,:,1,1,2]*D1cCG[:,:,1] +
    dXdxIC[:,:,1,2,2]*D2cCG[:,:,1] +
    dXdxIC[:,:,1,3,2].*D3cCGE[:,:,1] )
  @inbounds for iz=2:nz-1
    @views mul!(D1cCG[:,:,iz],CG.DS,cCG[:,:,iz])
    @views mul!(D2cCG[:,:,iz],cCG[:,:,iz],CG.DST)
    @views @. D3cCGE[:,:,iz] = 0.5*(D3cCG[:,:,iz-1] + D3cCG[:,:,iz]);
    @views @. F[:,:,iz,uPos] -= (dXdxIC[:,:,iz,1,1]*D1cCG[:,:,iz] +
      dXdxIC[:,:,iz,2,1]*D2cCG[:,:,iz] +
      dXdxIC[:,:,iz,3,1]*D3cCGE[:,:,iz] )
    @views @. F[:,:,iz,vPos] -= (dXdxIC[:,:,iz,1,2]*D1cCG[:,:,iz] +
      dXdxIC[:,:,iz,2,2]*D2cCG[:,:,iz] +
      dXdxIC[:,:,iz,3,2]*D3cCGE[:,:,iz] )
  end
end

function FGrad3RhoVec!(F,cCG,RhoCG,CG,Global,iF)
    (;  uPos,
        vPos,
        wPos) = Global.Model
  OP=CG.OrdPoly+1;
  NF=Global.Grid.NumFaces;
  nz=Global.Grid.nz;

  D1cCG = Global.Cache.CacheC1
  D2cCG = Global.Cache.CacheC2
  D3cCG = Global.Cache.CacheC3
  D3cCGE = Global.Cache.CacheC4
  @views dXdxIF = Global.Metric.dXdxIF[:,:,:,:,:,iF];
  @views dXdxIC = Global.Metric.dXdxIC[:,:,:,:,:,iF];


  @inbounds for iz=1:nz-1
    @views @. D3cCG[:,:,iz] = 0.5*(cCG[:,:,iz+1] - cCG[:,:,iz])
    @views @. F[:,:,iz,wPos] -= dXdxIF[:,:,iz+1,3,3]*D3cCG[:,:,iz] / 
      (0.5*(RhoCG[:,:,iz]+RhoCG[:,:,iz+1]))
  end  
  @views mul!(D1cCG[:,:,1],CG.DS,cCG[:,:,1])
  @views mul!(D2cCG[:,:,1],cCG[:,:,1],CG.DST)
  if nz>1
    @views  @. D3cCGE[:,:,1] = 1.5*D3cCG[:,:,1] - 0.5*D3cCG[:,:,2]
    @views mul!(D1cCG[:,:,nz],CG.DS,cCG[:,:,nz])
    @views mul!(D2cCG[:,:,nz],cCG[:,:,nz],CG.DST)
    @views  @. D3cCGE[:,:,nz] = D3cCG[:,:,nz-1]
    @views @. F[:,:,nz,uPos] -= (dXdxIC[:,:,nz,1,1]*D1cCG[:,:,nz] +
      dXdxIC[:,:,nz,2,1]*D2cCG[:,:,nz] +
      dXdxIC[:,:,nz,3,1]*D3cCGE[:,:,nz] ) / RhoCG[:,:,nz]
    @views @. F[:,:,nz,vPos] -= (dXdxIC[:,:,nz,1,2]*D1cCG[:,:,nz] +
      dXdxIC[:,:,nz,2,2]*D2cCG[:,:,nz] +
      dXdxIC[:,:,nz,3,2]*D3cCGE[:,:,nz] ) / RhoCG[:,:,nz]
  else
    @views  @. D3cCGE[:,:,1] = 0.0  
  end  
  @views @. F[:,:,1,uPos] -= (dXdxIC[:,:,1,1,1]*D1cCG[:,:,1] +
    dXdxIC[:,:,1,2,1]*D2cCG[:,:,1] +
    dXdxIC[:,:,1,3,1].*D3cCGE[:,:,1] ) / RhoCG[:,:,1]
  @views @. F[:,:,1,vPos] -= (dXdxIC[:,:,1,1,2]*D1cCG[:,:,1] +
    dXdxIC[:,:,1,2,2]*D2cCG[:,:,1] +
    dXdxIC[:,:,1,3,2].*D3cCGE[:,:,1] ) / RhoCG[:,:,1]
  @inbounds for iz=2:nz-1
    @views mul!(D1cCG[:,:,iz],CG.DS,cCG[:,:,iz])
    @views mul!(D2cCG[:,:,iz],cCG[:,:,iz],CG.DST)
    @views @. D3cCGE[:,:,iz] = 0.5*(D3cCG[:,:,iz-1] + D3cCG[:,:,iz]);
    @views @. F[:,:,iz,uPos] -= (dXdxIC[:,:,iz,1,1]*D1cCG[:,:,iz] +
      dXdxIC[:,:,iz,2,1]*D2cCG[:,:,iz] +
      dXdxIC[:,:,iz,3,1].*D3cCGE[:,:,iz] ) / RhoCG[:,:,iz]
    @views @. F[:,:,iz,vPos] -= (dXdxIC[:,:,iz,1,2]*D1cCG[:,:,iz] +
      dXdxIC[:,:,iz,2,2]*D2cCG[:,:,iz] +
      dXdxIC[:,:,iz,3,2]*D3cCGE[:,:,iz] ) / RhoCG[:,:,iz]
  end
end


function FGrad3Vec1!(F,cCG,CG,Global)
    (;  uPos,
        vPos,
        wPos) = Global.Model
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

D1cCG = Global.Cache.CacheC1
D2cCG = Global.Cache.CacheC2
D3cCG = Global.Cache.CacheC3
D3cCGE = Global.Cache.CacheC4
dXdxIF = Global.Metric.dXdxIF
dXdxIC = Global.Metric.dXdxIC


@inbounds for iF=1:NF
  @inbounds for iz=1:nz-1
    @views @. D3cCG[:,:,iz,iF] = 0.5*(cCG[:,:,iz+1,iF] - cCG[:,:,iz,iF])
    @views @. F[:,:,iz,iF,wPos] -= dXdxIF[:,:,iz+1,3,3,iF]*D3cCG[:,:,iz,iF] 
  end  
  @views mul!(D1cCG[:,:,1,iF],CG.DS,cCG[:,:,1,iF])
  @views mul!(D2cCG[:,:,1,iF],cCG[:,:,1,iF],CG.DST)
  @views  @. D3cCGE[:,:,1,iF] = 1.5*D3cCG[:,:,1,iF] - 0.5*D3cCG[:,:,iF,2]
  @views @. F[:,:,1,iF,uPos] -= (dXdxIC[:,:,1,1,1,iF]*D1cCG[:,:,1,iF] +
    dXdxIC[:,:,1,2,1,iF]*D2cCG[:,:,1,iF] +
    dXdxIC[:,:,1,3,1,iF].*D3cCGE[:,:,1,iF] ) 
  @views @. F[:,:,1,iF,vPos] -= (dXdxIC[:,:,1,1,2,iF]*D1cCG[:,:,1,iF] +
    dXdxIC[:,:,1,2,2,iF]*D2cCG[:,:,1,iF] +
    dXdxIC[:,:,1,3,2,iF].*D3cCGE[:,:,1,iF] )
  @inbounds for iz=2:nz-1
    @views mul!(D1cCG[:,:,iz,iF],CG.DS,cCG[:,:,iz,iF])
    @views mul!(D2cCG[:,:,iz,iF],cCG[:,:,iz,iF],CG.DST)
    @views @. D3cCGE[:,:,:,iz] .= 0.5.*(D3cCG[:,:,:,iz-1] .+ D3cCG[:,:,:,iz]);
    @views @. F[:,:,iz,iF,uPos] -= (dXdxIC[:,:,iz,1,1,iF]*D1cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,2,1,iF]*D2cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,3,1,iF].*D3cCGE[:,:,iz,iF] )
    @views @. F[:,:,iz,iF,vPos] -= (dXdxIC[:,:,iz,1,2,iF]*D1cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,2,2,iF]*D2cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,3,2,iF].*D3cCGE[:,:,iz,iF] )
  end
  @views mul!(D1cCG[:,:,nz,iF],CG.DS,cCG[:,:,nz,iF])
  @views mul!(D2cCG[:,:,nz,iF],cCG[:,:,nz,iF],CG.DST)
  @views  @. D3cCGE[:,:,nz,iF] = D3cCG[:,:,nz-1,iF]
  @views @. F[:,:,nz,iF,uPos] -= (dXdxIC[:,:,nz,1,1,iF]*D1cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,2,1,iF]*D2cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,3,1,iF].*D3cCGE[:,:,nz,iF] )
  @views @. F[:,:,nz,iF,vPos] -= (dXdxIC[:,:,nz,1,2,iF]*D1cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,2,2,iF]*D2cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,3,2,iF].*D3cCGE[:,:,nz,iF] )
end  
end

function FGrad3RhoVec1!(F,cCG,RhoCG,CG,Global)
    (;  uPos,
        vPos,
        wPos) = Global.Model
OP=CG.OrdPoly+1;
NF=Global.Grid.NumFaces;
nz=Global.Grid.nz;

D1cCG = Global.Cache.CacheC1
D2cCG = Global.Cache.CacheC2
D3cCG = Global.Cache.CacheC3
D3cCGE = Global.Cache.CacheC4
dXdxIF = Global.Metric.dXdxIF
dXdxIC = Global.Metric.dXdxIC


@inbounds for iF=1:NF
  for iz=1:nz-1
    @views @. D3cCG[:,:,iz,iF] = 0.5*(cCG[:,:,iz+1,iF] - cCG[:,:,iz,iF])
    @views @. F[:,:,iz,iF,wPos] -= dXdxIF[:,:,iz+1,3,3,iF]*D3cCG[:,:,iz,iF] / 
      (0.5*(RhoCG[:,:,iz,iF]+RhoCG[:,:,iz+1,iF]))
  end  
  @views mul!(D1cCG[:,:,1,iF],CG.DS,cCG[:,:,1,iF])
  @views mul!(D2cCG[:,:,1,iF],cCG[:,:,1,iF],CG.DST)
  @views  @. D3cCGE[:,:,1,iF] = 1.5*D3cCG[:,:,1,iF] - 0.5*D3cCG[:,:,iF,2]
  @views @. F[:,:,1,iF,uPos] -= (dXdxIC[:,:,1,1,1,iF]*D1cCG[:,:,1,iF] +
    dXdxIC[:,:,1,2,1,iF]*D2cCG[:,:,1,iF] +
    dXdxIC[:,:,1,3,1,iF].*D3cCGE[:,:,1,iF] ) / RhoCG[:,:,1,iF]
  @views @. F[:,:,1,iF,vPos] -= (dXdxIC[:,:,1,1,2,iF]*D1cCG[:,:,1,iF] +
    dXdxIC[:,:,1,2,2,iF]*D2cCG[:,:,1,iF] +
    dXdxIC[:,:,1,3,2,iF].*D3cCGE[:,:,1,iF] ) / RhoCG[:,:,1,iF]
  @inbounds for iz=2:nz-1
    @views mul!(D1cCG[:,:,iz,iF],CG.DS,cCG[:,:,iz,iF])
    @views mul!(D2cCG[:,:,iz,iF],cCG[:,:,iz,iF],CG.DST)
    @views @. D3cCGE[:,:,:,iz] .= 0.5.*(D3cCG[:,:,:,iz-1] .+ D3cCG[:,:,:,iz]);
    @views @. F[:,:,iz,iF,uPos] -= (dXdxIC[:,:,iz,1,1,iF]*D1cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,2,1,iF]*D2cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,3,1,iF].*D3cCGE[:,:,iz,iF] ) / RhoCG[:,:,iz,iF]
    @views @. F[:,:,iz,iF,vPos] -= (dXdxIC[:,:,iz,1,2,iF]*D1cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,2,2,iF]*D2cCG[:,:,iz,iF] +
      dXdxIC[:,:,iz,3,2,iF].*D3cCGE[:,:,iz,iF] ) / RhoCG[:,:,iz,iF]
  end
  @views mul!(D1cCG[:,:,nz,iF],CG.DS,cCG[:,:,nz,iF])
  @views mul!(D2cCG[:,:,nz,iF],cCG[:,:,nz,iF],CG.DST)
  @views  @. D3cCGE[:,:,nz,iF] = D3cCG[:,:,nz-1,iF]
  @views @. F[:,:,nz,iF,uPos] -= (dXdxIC[:,:,nz,1,1,iF]*D1cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,2,1,iF]*D2cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,3,1,iF].*D3cCGE[:,:,nz,iF] ) / RhoCG[:,:,nz,iF]
  @views @. F[:,:,nz,iF,vPos] -= (dXdxIC[:,:,nz,1,2,iF]*D1cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,2,2,iF]*D2cCG[:,:,nz,iF] +
    dXdxIC[:,:,nz,3,2,iF].*D3cCGE[:,:,nz,iF] ) / RhoCG[:,:,nz,iF]
end  
end


