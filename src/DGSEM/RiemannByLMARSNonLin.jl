@inline function RiemannByLMARSNonLin!(F,VLL,VRR)
  hPos = 1
  uPos = 2
  vPos = 3

  cS = sqrt(9.81 * 1.e5)

  pLL = PresSh(VLL)
  pRR = PresSh(VRR)
  hM = 0.5 * (VLL[hPos] + VRR[hPos])
  vLL = VLL[uPos] / VLL[hPos]
  vRR = VRR[uPos] / VRR[hPos]
  pM = 0.5 * (pLL + pRR) - 0.5 * cS * hM * (vRR - vLL)
  vM = 0.5 * (vRR + vLL) -1.0 /(2.0 * cS) * (pRR - pLL) / hM
  if vM > 0
    F[hPos] = vM * VLL[hPos]
    F[uPos] = vM * VLL[uPos] + pM
    F[vPos] = vM * VLL[vPos]
  else
    F[hPos] = vM * VRR[hPos]
    F[uPos] = vM * VRR[uPos] + pM
    F[vPos] = vM * VRR[vPos]
  end
end    
