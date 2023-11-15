function Limit!(qMin,qMax,Rhoq,Rho,CG,Global)
  OP=CG.OrdPoly+1
  NF=Global.Grid.NumFaces
  nz=Global.Grid.nz
  @. qMin = 1.e40
  @. qMax = -1.e40
  @inbounds for iT = 1:size(Rhoq,3)
    @inbounds for iF = 1:NF
      @inbounds for iD = 1:OP*OP
        ind = CG.Glob[iD,iF]
        @inbounds for iz=1:nz
          qMin[iz,iF,iT]=min(qMin[iz,iF,iT],Rhoq[iz,ind,iT]/Rho[iz,ind])
          qMax[iz,iF,iT]=max(qMax[iz,iF,iT],Rhoq[iz,ind,iT]/Rho[iz,ind])
        end
      end
    end
  end
end

function HorLimiter!(dcRho,cRho,RhoS,RhoSS,dt,J,w,qMin,qMax,ThreadCache,iF)
  @unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = ThreadCache
  l0 = 0
  eta = 1.e-12
  etaR = 1.e-10
  dlFD = 1.e-8
  n = size(w,1)
  WM = TCacheC1[Threads.threadid()]
  q = TCacheC2[Threads.threadid()]
  JC = TCacheC3[Threads.threadid()]
  cRhoS = TCacheC4[Threads.threadid()]
  SumJ = sum(J)
  for i = 1 : n
    for j = 1 : n
      iD = i + (j - 1) * n  
      JC[iD] = (J[iD,1] + J[iD,2])  
      WM[iD] = JC[iD] * w[i] * w[j] / SumJ
    end
  end
  @. cRhoS = cRho + dt * dcRho / JC
# Finite difference step
  resp = 0
  @inbounds for i in eachindex(q)
    q[i] = median([qMin,cRhoS[i]/RhoS[i]+l0*WM[i]*RhoSS[i]/RhoS[i],qMax])
    resp = resp + WM[i]*(q[i]*RhoSS[i]-cRhoS[i])
  end
  resp0 = abs(resp)
  if abs(resp) <= eta
    @. dcRho = (q * RhoSS - cRho) * JC / dt
    return
  end
  resc = 0
  @inbounds for i in eachindex(q)
    qLoc = median([qMin,cRhoS[i]/RhoS[i]+(l0+dlFD)*WM[i]*RhoSS[i]/RhoS[i],qMax])
    resc = resc + WM[i] * (qLoc * RhoSS[i] - cRhoS[i])
  end
  if abs(resc-resp) <= 1.e-13
    @. dcRho = (q * RhoSS - cRho) * JC / dt
    return
  end  
  alpha = dlFD / (resc - resp)
  lp = l0
  lc = lp - alpha * resp
  iter = 1
  resp = 0
  while abs(resc) > eta * resp0
    resc = 0
    @inbounds for i in eachindex(q)
      q[i] = median([qMin,cRhoS[i]/RhoS[i]+lc*WM[i]*RhoSS[i]/RhoS[i],qMax])
      resc = resc + WM[i]*(q[i]*RhoSS[i]-cRhoS[i])
    end
    if abs(resc-resp) <= 1.e-13
      @. dcRho = (q * RhoSS - cRho) * JC / dt
      return
    end  
    alpha = (lp - lc) / (resp - resc)
    resp = resc
    lp = lc
    lc = lc - alpha * resc
    iter = iter + 1
  end
  @. dcRho = (q * RhoSS - cRho) * JC / dt
end 

