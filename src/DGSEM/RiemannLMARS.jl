function (::RiemannLMARS)(Param,Phys,RhoPos,uPos,vPos,wPos,ThPos,pPos)
  @inline function RiemannByLMARSNonLin!(F,VLL,VRR,AuxL,AuxR,Normal)

    FT = eltype(F)
    cS = Param.cS
    pLL = AuxL[pPos]
    pRR = AuxR[pPos]
    RhoM = FT(0.5) * (VLL[RhoPos] + VRR[RhoPos])
    vLL = (VLL[uPos] * Normal[1] + VLL[vPos] * Normal[2] + VLL[wPos] * Normal[3]) / VLL[RhoPos]
    vRR = (VRR[uPos] * Normal[1] + VRR[vPos] * Normal[2] + VRR[wPos] * Normal[3]) / VRR[RhoPos]
    pM = FT(0.5) * (pLL + pRR) - FT(0.5) * cS * RhoM * (vRR - vLL)
    vM = FT(0.5) * (vRR + vLL) - FT(1.0) /(FT(2.0) * cS) * (pRR - pLL) / RhoM
    if vM > FT(0)
      F[RhoPos] = vM * VLL[RhoPos]
      F[uPos] = vM * VLL[uPos] + Normal[1] * pM
      F[vPos] = vM * VLL[vPos] + Normal[2] * pM
      F[wPos] = vM * VLL[wPos] + Normal[3] * pM
      F[ThPos] = vM * VLL[ThPos]
    else
      F[RhoPos] = vM * VRR[RhoPos]
      F[uPos] = vM * VRR[uPos] + Normal[1] * pM
      F[vPos] = vM * VRR[vPos] + Normal[2] * pM
      F[wPos] = vM * VRR[wPos] + Normal[3] * pM
      F[ThPos] = vM * VRR[ThPos]
    end
  end
  return RiemannByLMARSNonLin!
end
