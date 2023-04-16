function Limit!(qMin,qMax,Rhoq,Rho,CG,Global)
  OP=CG.OrdPoly+1
  NF=Global.Grid.NumFaces
  nz=Global.Grid.nz
  @. qMin = 1.e40
  @. qMax = -1.e40
  @inbounds for iT = 1:size(Rhoq,3)
    @inbounds for iF = 1:NF
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            qMin[iz,iF,iT]=min(qMin[iz,iF],Rhoq[iz,ind,iT]/Rho[iz,ind])
            qMax[iz,iF,iT]=max(qMax[iz,iF],Rhoq[iz,ind,iT]/Rho[iz,ind])
          end
        end
      end
    end
  end
end

function HorLimiter!(dcRho,cRho,RhoS,RhoSS,dt,J,w,qMin,qMax)
  l0 = 0
  eta = 1.e-12
  dlFD = 1.e-8
  n = size(w,1)
  WM = zeros(n,n)
  q = zeros(n,n)
  SumJ = sum(J)
  for i = 1 : n
    for j = 1 : n
      WM[i,j] = J[i,j] * w[i] * w[j] / SumJ
    end
  end
  cRhoS = zeros(n,n)
  @. cRhoS = cRho + dt * dcRho / J
# Finite difference step
  resp = 0
  @inbounds for i in eachindex(q)
    q[i] = median([qMin,cRhoS[i]/RhoS[i]+l0*WM[i]*RhoSS[i]/RhoS[i],qMax])
    resp = resp + WM[i]*(q[i]*RhoSS[i]-cRhoS[i])
  end
  if abs(resp) <= eta
    @. dcRho = (q * RhoSS - cRho) * J / dt
    return
  end
  resc = 0
  @inbounds for i in eachindex(q)
    qLoc = median([qMin,cRhoS[i]/RhoS[i]+(l0+dlFD)*WM[i]*RhoSS[i]/RhoS[i],qMax])
    resc = resc + WM[i] * (qLoc * RhoSS[i] - cRhoS[i])
  end
  if abs(resc-resp) <= 1.e-13
    @. dcRho = (q * RhoSS - cRho) * J / dt
    return
  end  
  alpha = dlFD / (resc - resp)
  lp = l0
  lc = lp - alpha * resp
  iter = 0
  while abs(resc) > eta
    resc = 0
    @inbounds for i in eachindex(q)
      q[i] = median([qMin,cRhoS[i]/RhoS[i]+lc*WM[i]*RhoSS[i]/RhoS[i],qMax])
      resc = resc + WM[i]*(q[i]*RhoSS[i]-cRhoS[i])
    end
    if abs(resc-resp) <= 1.e-13
      @. dcRho = (q * RhoSS - cRho) * J / dt
      return
    end  
    alpha = (lp - lc) / (resp - resc)
    resp = resc
    lp = lc
    lc = lc - alpha * resc
    iter = iter + 1
  end
  @. dcRho = (q * RhoSS - cRho) * J / dt
end 

function QP!(dcRho,cRho,RhoS,RhoSS,dt,J,w,qMin,qMax)
# QStar (m^3/kg)
# Rho  (m^3/kg) 
  n = size(w,1)
  WM = zeros(n,n)
  cRhoS = zeros(n,n)
  q = zeros(n,n)
  JRho = zeros(n,n)

  for i = 1 : n
    for j = 1 : n
      WM[i,j] = J[i,j] * w[i] * w[j] 
    end
  end
  @. cRhoS = cRho + dt * dcRho / J

  tol_limiter=5e-14
  MassRho=0.0
  Massq=0.0
  for i in eachindex(cRhoS) 
    JRho[i] = J[i] * RhoS[i]
    q[i] = cRhoS[i] / RhoS[i]
    MassRho = MassRho + JRho[i] 
    Massq = Massq + JRho[i] * q[i] 
  end

  if MassRho <= 0
    @show "MassRho <= 0"  
    return
  end
  qMinOld = qMin
  qMaxOld = qMax

# relax constraints to ensure limiter has a solution:
# This is only needed if running with the SSP CFL>1 or
# due to roundoff errors
  if Massq < qMin*MassRho
    qMin=Massq/MassRho
  end
  if Massq > qMax*MassRho
    qMax=Massq/MassRho
  end
  for iter=1:10
    addMassq=0.0
#   Computation of the error tolerance 
#   and projecting x into the upper and
#   lower bounds
    for i in eachindex(cRhoS)
      if q[i]>qMax
        addMassq=addMassq+(q[i]-qMax)*JRho[i]
        q[i]=qMax
      elseif q[i]<qMin
        addMassq=addMassq-(qMin-q[i])*JRho[i]
        q[i]=qMin
      end
    end
  
    weightssum=0.0
    if addMassq>0
      for i in eachindex(cRhoS)
        if q[i]<qMax
          weightssum=weightssum+JRho[i]
        end
      end
      for i in eachindex(cRhoS)
        if q[i]<qMax
          q[i]=q[i]+addMassq/weightssum
        end
      end
    else
      for i in eachindex(cRhoS)
        if q[i]>qMin
          weightssum=weightssum+JRho[i]
        end
      end
      for i in eachindex(cRhoS)
        if q[i]>qMin
          q[i]=q[i]+addMassq/weightssum
        end
      end
    end
  end
# @. cRhoS = cRho + dt * dcRho / J
  @. dcRho = (q*RhoS - cRho) * J / dt
end   

