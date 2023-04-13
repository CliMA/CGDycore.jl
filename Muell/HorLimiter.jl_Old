mutable struct Stencil
  F::Array{Int,1}
  P::Array{Int,1}
end

function Stencil()
  F = zeros(0)
  P = zeros(0)
  return Stencil(
    F,
    P
  )
end  

function Stencil(Grid::GridStruct,iF)
  FLoc=zeros(Int,16)
  PLoc=zeros(Int,16)
  iS=0
  for i=1:4
    iN=Grid.Faces[iF].N[i]
    for j=1:size(Grid.Nodes[iN].FG,1)
      jF=Grid.Nodes[iN].FG[j]
      jP=Grid.Nodes[iN].FP[j]
      inside=false
      for jS=1:iS
        if FLoc[jS]==jF
          inside=true
          break
        end
      end
      if !inside
        iS=iS+1
        FLoc[iS]=jF
        PLoc[iS]=jP
      end
    end
  end
  F = zeros(Int,iS)
  P = zeros(Int,iS)
  F .= FLoc[1:iS]
  P .= PLoc[1:iS]
  return Stencil(
    F,
    P
  )
end  

mutable struct HorLimiter
  NumStencils::Int
  Stencils::Array{Stencil, 1}
end

function HorLimiter(Grid)
  NumStencils = Grid.NumFaces
  Stencils = map(1:NumStencils) do i
    Stencil(Grid,i)
  end  
  return HorLimiter(
    NumStencils,
    Stencils
  )
end

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

function HorLimiter!(RhoqNew,Fac1,Rhoq,Fac2,dRhoq,Rho,CG,Global)
  OP=CG.OrdPoly+1
  NF=Global.Grid.NumFaces
  nz=Global.Grid.nz
  cLoc = Global.Cache.uStar
  @views qMin = Global.Cache.Pres[1,1,:,:]
  @views qMax = Global.Cache.Pres[1,2,:,:]
  qStarCG = Global.Cache.CacheE1
  RhoCG = Global.Cache.CacheE2
  cCorrCG = Global.Cache.CacheE3
  JCMass = Global.Cache.CacheE4
  Tr = zeros(OP,OP)
  dTr = zeros(OP,OP)
  @inbounds for iT = 1:size(RhoqNew,3)
    @. qMin = 1.e40
    @. qMax = -1.e40
    @inbounds for iF = 1:NF
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            qMin[iz,iF]=min(qMin[iz,iF],Rhoq[iz,ind,iT]/Rho[iz,ind])
            qMax[iz,iF]=max(qMax[iz,iF],Rhoq[iz,ind,iT]/Rho[iz,ind])
          end
        end
      end
    end
    @inbounds for iF=1:NF
      @inbounds for iz=1:nz  
        @views qMinS=minimum(qMin[iz,CG.Stencil[iF,:]])
        @views qMaxS=maximum(qMax[iz,CG.Stencil[iF,:]])
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
            qStarCG[iP,jP]=Rhoq[iz,ind,iT]+Fac2*dRhoq[iz,ind,iT]
            RhoCG[iP,jP] = Rho[iz,ind]
            JCMass[iP,jP] = Global.Metric.JC[iP,jP,iz,iF] * CG.w[iP] * CG.w[jP]
          end
        end  
        QP!(cCorrCG,qStarCG,JCMass,qMinS,qMaxS,RhoCG)
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
            RhoqNew[iz,ind,iT] += Fac1*cCorrCG[iP,jP] / CG.MMass[iz,ind]
          end
        end
      end
    end
  end
end

function QP!(RhoQCorr,RhoQStar,J,QMin,QMax,Rho)
# QStar (m^3/kg)
# Rho  (m^3/kg) 
  tol_limiter=5e-14
  JRho = J
  Q = RhoQCorr
  MassRho=0.0
  MassQ=0.0
  @inbounds for i in eachindex(RhoQStar) 
    JRho[i]=J[i]*Rho[i]
    Q[i]=RhoQStar[i]/Rho[i] # Mixing ratio
    MassRho=MassRho+JRho[i] # Mass of density
    MassQ=MassQ+JRho[i]*Q[i] # Mass of the tracer 
  end

  if MassRho <= 0
    return
  end

  QMinOld=QMin
  QMaxOld=QMax


# relax constraints to ensure limiter has a solution:
# This is only needed if running with the SSP CFL>1 or
# due to roundoff errors
  if MassQ < QMin*MassRho
    QMin=MassQ/MassRho
  end
  if MassQ > QMax*MassRho
    QMax=MassQ/MassRho
  end
  @inbounds for iter=1:5
    addMassQ=0.0
#   Computation of the error tolerance 
#   and projecting x into the upper and
#   lower bounds
    @inbounds for i in eachindex(RhoQStar)
      if Q[i]>QMax
        addMassQ=addMassQ+(Q[i]-QMax)*JRho[i]
        Q[i]=QMax
      elseif Q[i]<QMin
        addMassQ=addMassQ-(QMin-Q[i])*JRho[i]
        Q[i]=QMin
      end
    end
  
    weightssum=0.0
    if addMassQ>0
      @inbounds for i in eachindex(RhoQStar)
        if Q[i]<QMax
          weightssum=weightssum+JRho[i]
        end
      end
      @inbounds for i in eachindex(RhoQStar)
        if Q[i]<QMax
          Q[i]=Q[i]+addMassQ/weightssum
        end
      end
    else
      @inbounds for i in eachindex(RhoQStar)
        if Q[i]>QMin
          weightssum=weightssum+JRho[i]
        end
      end
      @inbounds for i in eachindex(RhoQStar)
        if Q[i]>QMin
          Q[i]=Q[i]+addMassQ/weightssum
        end
      end
    end
  end

  @. RhoQCorr = Q * JRho
end   

function QPdt!(dTr,Tr,dRho,Rho,QMin,QMax,dt,J,w)
# QStar (m^3/kg)
# Rho  (m^3/kg)
  tol_limiter=5e-14
  OP = size(dTr,1)
  RhoNew = zeros(OP,OP)
  JRho = zeros(OP,OP)
  MassRho=0.0
  MassQ=0.0
  Q = dTr 
  for j = 1 : OP
    for i = 1 : OP
      RhoNew[i,j] = Rho[i,j] + dt * dRho[i,j] / J[i,j]  
      JRho[i,j] = RhoNew[i,j] * J[i,j] * w[i] * w[j]  
      Q[i,j] = (Tr[i,j] * J[i,j] + dt * dTr[i,j]) / RhoNew[i,j] / J[i,j]
      MassRho += JRho[i,j]
      MassQ += Q[i,j] * JRho[i,j]
    end
  end  

  QMinOld=QMin
  QMaxOld=QMax

  # relax constraints to ensure limiter has a solution:
# This is only needed if running with the SSP CFL>1 or
# due to roundoff errors
  if MassQ < QMin*MassRho
    QMin=MassQ/MassRho
  end
  if MassQ > QMax*MassRho
    QMax=MassQ/MassRho
  end
  @inbounds for iter=1:5
    addMassQ=0.0
    #   Computation of the error tolerance
#   and projecting x into the upper and
#   lower bounds
    @inbounds for i in eachindex(Q)
      if Q[i]>QMax
        addMassQ=addMassQ+(Q[i]-QMax)*JRho[i]
        Q[i]=QMax
      elseif Q[i]<QMin
        addMassQ=addMassQ-(QMin-Q[i])*JRho[i]
        Q[i]=QMin
      end
    end
    weightssum=0.0
    if addMassQ>0
      @inbounds for i in eachindex(Q)
        if Q[i]<QMax
          weightssum=weightssum+JRho[i]
        end
      end
      @inbounds for i in eachindex(Q)
        if Q[i]<QMax
          Q[i]=Q[i]+addMassQ/weightssum
        end
      end
    else
      @inbounds for i in eachindex(Q)
        if Q[i]>QMin
          weightssum=weightssum+JRho[i]
        end
      end
      @inbounds for i in eachindex(Q)
        if Q[i]>QMin
          Q[i]=Q[i]+addMassQ/weightssum
        end
      end
    end
  end
  @. dTr = (Q * RhoNew - Tr) * J / dt 
# Q * RhoNew = Tr + dt * dTr / J 

end

function QPSecant!(dTr,Tr,dRho,Rho,QMin,QMax,dt,J,w)
# QStar (m^3/kg)
# Rho  (m^3/kg)
  tol_limiter=5e-14
  OP = size(dTr,1)
  RhoNew = zeros(OP,OP)
  w = RhoNew
  JRho = zeros(OP,OP)
  W = JRho
  dm = zeros(OP,OP)
  Q = dTr
  eta = 1.e-12
  dlFD = 1.e-6
  @inbounds for j = 1 : OP
    @inbounds for i = 1 : OP
      RhoNew[i,j] = Rho[i,j] + dt * dRho[i,j] / J[i,j]
      JRho[i,j] = RhoNew[i,j] * J[i,j] * w[i] * w[j]
      Q[i,j] = (Tr[i,j] * J[i,j] + dt * dTr[i,j]) / RhoNew[i,j] / J[i,j]
    end
  end
  @show RhoNew
  @show Q
  @show QMin,QMax

  l0 = 0.0
  resp = 0.0
  @inbounds for i in eachindex(Q)
    dm[i] = median([QMin - Q[i],l0*W[i]/w[i],QMax - Q[i]])
    resp = resp + W[i] * dm[i]
  end
  if abs(resp) <= eta
   return
  end
  resc = 0.0
  @inbounds for i in eachindex(Q)
    dm[i] = median([QMin - Q[i],(l0+dlFD)*W[i]/w[i],QMax - Q[i]]);
    resc = resc + W[i] * dm[i]
  end
  alpha = dlFD / (resc - resp) 
  lp = l0 
  lc = lp - alpha * resp 
  @show "Start",lc
  while abs(resc) > eta
    resc = 0
    @inbounds for i in eachindex(Q)
      dm[i] = median([QMin - Q[i],lc*W[i]/w[i],QMax - Q[i]]);
      resc = resc + W[i] * dm[i]
    end
    @show resp,resc
    @show lp,lc
    alpha = (lp - lc)/(resp - resc)
    resp = resc
    lp = lc
    lc = lc - alpha * resc
    @show lc
    @show dm
  end
  @. Q = Q + dm
  @. dTr = (Q * RhoNew - Tr) * J / dt 
end  

#function dm = SecantWMethod(w,W,dmMin,dmMax)
#l0 = 0;
#eta = 1.e-12;
#dlFD = 1.e-6;
#n = size(w,1);
#dm = zeros(size(w));
#% Finite difference step
#resp = 0;
#for i = 1 : n
#  dm(i) = median([dmMin(i),l0*W(i)/w(i),dmMax(i)]);
#  resp = resp + W(i)*dm(i);
#end
#if abs(resp) <= eta
#  return
#end
#resc = 0;
#for i = 1 : n
#  dm(i) = median([dmMin(i),(l0+dlFD)*W(i)/w(i),dmMax(i)]);
#  resc = resc + W(i)*dm(i);
#end
#alpha = dlFD / (resc - resp);
#lp = l0;
#lc = lp - alpha * resp;
#while abs(resc) > eta
#  resc = 0;
#  for i = 1 :n
#    dm(i) = median([dmMin(i),lc*W(i)/w(i),dmMax(i)]);
#    resc = resc + W(i)*dm(i);
#  end
#  alpha = (lp - lc)/(resp - resc);
#  resp = resc;
#  lp = lc;
#  lc = lc - alpha * resc
#end
#
#end
  

function SecantQ!(dcRho,cRho,RhoS,RhoSS,dt,J,w,qMin,qMax)
  l0 = 0
  eta = 1.e-10
  dlFD = 1.e-8
  n = size(w,1)
  WM = zeros(n,n)
  q = zeros(n,n)
  for i = 1 : n
    for j = 1 : n
      WM[i,j] = J[i,j] * w[i] * w[j]
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
    alpha = (lp - lc) / (resp - resc)
    resp = resc
    lp = lc
    lc = lc - alpha * resc
    iter = iter + 1
#   @show iter,resc
  end
  @. dcRho = (q * RhoSS - cRho) * J / dt
end 
