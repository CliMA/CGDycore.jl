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
  Q = Rho
  MassRho=0.0
  MassQ=0.0
  for i in eachindex(RhoQStar) 
    JRho[i]=J[i]*Rho[i]
    Q[i]=RhoQStar[i]/Rho[i]
    MassRho=MassRho+JRho[i] 
    MassQ=MassQ+JRho[i]*Q[i] 
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
  for iter=1:5
    addMassQ=0.0
#   Computation of the error tolerance 
#   and projecting x into the upper and
#   lower bounds
    for i in eachindex(RhoQStar)
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
      for i in eachindex(RhoQStar)
        if Q[i]<QMax
          weightssum=weightssum+JRho[i]
        end
      end
      for i in eachindex(RhoQStar)
        if Q[i]<QMax
          Q[i]=Q[i]+addMassQ/weightssum
        end
      end
    else
      for i in eachindex(RhoQStar)
        if Q[i]>QMin
          weightssum=weightssum+JRho[i]
        end
      end
      for i in eachindex(RhoQStar)
        if Q[i]>QMin
          Q[i]=Q[i]+addMassQ/weightssum
        end
      end
    end
  end
  RhoQCorr .= Q.*JRho
end   
