function (::KennedyGruberGrav)(RhoPos,uPos,vPos,wPos,ThPos,pPos,GPPos)
  @inline function FluxNonLinAver!(flux,VL,VR,AuxL,AuxR,m_L,m_R)
    FT = eltype(flux)

    pAv = FT(0.5) * ((AuxL[pPos] + AuxR[pPos]) + 
      FT(0.5) * (VL[RhoPos] + VR[RhoPos]) * (AuxR[GPPos] - AuxL[GPPos]))
    uAv = FT(0.5) * (VL[uPos] + VR[uPos])
    vAv = FT(0.5) * (VL[vPos] + VR[vPos])
    wAv = FT(0.5) * (VL[wPos] + VR[wPos])
    RhoAv = FT(0.5) * (VL[RhoPos] + VR[RhoPos])
    ThAv = FT(0.5) * (VL[ThPos] + VR[ThPos])
    mAv1 = FT(0.5) * (m_L[1] + m_R[1])
    mAv2 = FT(0.5) * (m_L[2] + m_R[2])
    mAv3 = FT(0.5) * (m_L[3] + m_R[3])
    qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
    flux[1] = RhoAv * qHat
    flux[2] = flux[1] * uAv + mAv1 * pAv
    flux[3] = flux[1] * vAv + mAv2 * pAv
    flux[4] = flux[1] * wAv + mAv3 * pAv
    flux[5] = flux[1] * ThAv

  end
  return FluxNonLinAver!
end
