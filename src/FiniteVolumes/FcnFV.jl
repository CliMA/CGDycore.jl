function CacheFV(backend,FTB,MetricFV,Grid)
  CurlUu = zeros(FTB,Grid.NumNodes)
  TangUu = zeros(FTB,Grid.NumEdges)
  hE = zeros(FTB,Grid.NumEdges)
  K = zeros(FTB,Grid.NumFaces)
  Grad,Inter = GradMPFATri(backend,FTB,Grid)
  Div = DivMPFA(backend,FTB,MetricFV,Grid)
  Curl = CurlNodeMatrix(MetricFV,Grid)
  Tang = TagentialVelocity2Matrix(MetricFV,Grid)
  return CacheFV(
  CurlUu,
  TangUu,
  hE,
  K,
  Grad,
  Inter,
  Div,
  Curl,
  Tang,
    )
end  

function FcnFV!(F,U,MetricFV,Grid,Cache,Phys)

  pPosS = 1
  pPosEI = Grid.NumFaces
  pPosE = Grid.NumFaces + Grid.NumFacesG
  uPosS = pPosE + 1
  uPosE = pPosE + Grid.NumEdges
  @views rp = F[pPosS:pPosE]
  @views rpI = F[pPosS:pPosEI]
  @views ru = F[uPosS:uPosE]
  @views Up = U[pPosS:pPosE]
  @views UpI = U[pPosS:pPosEI]
  @views Uu = U[uPosS:uPosE]

  CurlUu = Cache.CurlUu
  @views TangUu = Cache.TangUu
  hE = Cache.hE
  grad = Cache.hE
  K = Cache.K
  Curl = Cache.Curl
  Tang = Cache.Tang
  Grad = Cache.Grad
  Div = Cache.Div
  Inter = Cache.Inter

  mul!(TangUu,Tang,Uu)
  KineticEnergy(K,Uu,TangUu,MetricFV,Grid)
  @. K += Up * Phys.Grav
  mul!(CurlUu,Curl,Uu)
  mul!(grad,Grad,K)
  @inbounds for iE = 1 : Grid.NumEdges
    iN1 = Grid.Edges[iE].N[1]
    iN2 = Grid.Edges[iE].N[2]
    x = Grid.Edges[iE].Mid.x
    y = Grid.Edges[iE].Mid.y
    z = Grid.Edges[iE].Mid.z
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    CurlEdge = (CurlUu[iN1] * MetricFV.DualVolume[iN1] + CurlUu[iN2] * MetricFV.DualVolume[iN2]) / 
      (MetricFV.DualVolume[iN1] + MetricFV.DualVolume[iN2])  
    ru[iE] = -TangUu[iE] * (CurlEdge + 2 * Phys.Omega * sinlat) - grad[iE]
  end
  mul!(hE,Inter,Up)
  @. hE *= Uu
  mul!(rpI,Div,hE)
end
function FcnFV1!(F,U,MetricFV,Grid,Cache,Phys)

  pPosS = 1
  pPosEI = Grid.NumFaces
  pPosE = Grid.NumFaces + Grid.NumFacesG
  uPosS = pPosE + 1
  uPosE = pPosE + Grid.NumEdges
  @views rp = F[pPosS:pPosE]
  @views rpI = F[pPosS:pPosEI]
  @views ru = F[uPosS:uPosE]
  @views Up = U[pPosS:pPosE]
  @views UpI = U[pPosS:pPosEI]
  @views Uu = U[uPosS:uPosE]

  CurlUu = Cache.CurlUu
  @views TangUu = Cache.TangUu
  hE = Cache.hE
  grad = Cache.hE
  K = Cache.K
  Curl = Cache.Curl
  Tang = Cache.Tang
  Grad = Cache.Grad
  Div = Cache.Div
  Inter = Cache.Inter

  mul!(TangUu,Tang,Uu)
  KineticEnergy(K,Uu,TangUu,MetricFV,Grid)
  @. K += Up * Phys.Grav
  mul!(CurlUu,Curl,Uu)
  mul!(grad,Grad,K)
  @inbounds for iE = 1 : Grid.NumEdges
    iN1 = Grid.Edges[iE].N[1]
    iN2 = Grid.Edges[iE].N[2]
    x = Grid.Edges[iE].Mid.x
    y = Grid.Edges[iE].Mid.y
    z = Grid.Edges[iE].Mid.z
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    CurlEdge = (CurlUu[iN1] * MetricFV.DualVolume[iN1] + CurlUu[iN2] * MetricFV.DualVolume[iN2]) / 
      (MetricFV.DualVolume[iN1] + MetricFV.DualVolume[iN2])  
    ru[iE] = -TangUu[iE] * (CurlEdge + 2 * Phys.Omega * sinlat) - grad[iE]
    ru[iE] = - grad[iE]
  end
  mul!(hE,Inter,Up)
  @. hE *= Uu
  mul!(rpI,Div,hE)
end
function FcnFV2!(F,U,MetricFV,Grid,Cache,Phys)

  pPosS = 1
  pPosEI = Grid.NumFaces
  pPosE = Grid.NumFaces + Grid.NumFacesG
  uPosS = pPosE + 1
  uPosE = pPosE + Grid.NumEdges
  @views rp = F[pPosS:pPosE]
  @views rpI = F[pPosS:pPosEI]
  @views ru = F[uPosS:uPosE]
  @views Up = U[pPosS:pPosE]
  @views UpI = U[pPosS:pPosEI]
  @views Uu = U[uPosS:uPosE]

  CurlUu = Cache.CurlUu
  @views TangUu = Cache.TangUu
  hE = Cache.hE
  grad = Cache.hE
  K = Cache.K
  Curl = Cache.Curl
  Tang = Cache.Tang
  Grad = Cache.Grad
  Div = Cache.Div
  Inter = Cache.Inter

  mul!(TangUu,Tang,Uu)
  KineticEnergy(K,Uu,TangUu,MetricFV,Grid)
  @. K += Up * Phys.Grav
  mul!(CurlUu,Curl,Uu)
  mul!(grad,Grad,K)
  @inbounds for iE = 1 : Grid.NumEdges
    iN1 = Grid.Edges[iE].N[1]
    iN2 = Grid.Edges[iE].N[2]
    x = Grid.Edges[iE].Mid.x
    y = Grid.Edges[iE].Mid.y
    z = Grid.Edges[iE].Mid.z
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    CurlEdge = (CurlUu[iN1] * MetricFV.DualVolume[iN1] + CurlUu[iN2] * MetricFV.DualVolume[iN2]) / 
      (MetricFV.DualVolume[iN1] + MetricFV.DualVolume[iN2])  
    ru[iE] = -TangUu[iE] * (CurlEdge + 2 * Phys.Omega * sinlat) - grad[iE]
    ru[iE] = +TangUu[iE] * (CurlEdge + 2 * Phys.Omega * sinlat)
  end
  mul!(hE,Inter,Up)
  @. hE *= Uu
  mul!(rpI,Div,hE)
end
function FcnFV3!(F,U,MetricFV,Grid,Cache,Phys)

  pPosS = 1
  pPosEI = Grid.NumFaces
  pPosE = Grid.NumFaces + Grid.NumFacesG
  uPosS = pPosE + 1
  uPosE = pPosE + Grid.NumEdges
  @views rp = F[pPosS:pPosE]
  @views rpI = F[pPosS:pPosEI]
  @views ru = F[uPosS:uPosE]
  @views Up = U[pPosS:pPosE]
  @views UpI = U[pPosS:pPosEI]
  @views Uu = U[uPosS:uPosE]

  CurlUu = Cache.CurlUu
  @views TangUu = Cache.TangUu
  hE = Cache.hE
  grad = Cache.hE
  K = Cache.K
  Curl = Cache.Curl
  Tang = Cache.Tang
  Grad = Cache.Grad
  Div = Cache.Div
  Inter = Cache.Inter

  mul!(TangUu,Tang,Uu)
  KineticEnergy(K,Uu,TangUu,MetricFV,Grid)
  @. K += Up * Phys.Grav
  mul!(CurlUu,Curl,Uu)
  mul!(grad,Grad,K)
  @inbounds for iE = 1 : Grid.NumEdges
    iN1 = Grid.Edges[iE].N[1]
    iN2 = Grid.Edges[iE].N[2]
    x = Grid.Edges[iE].Mid.x
    y = Grid.Edges[iE].Mid.y
    z = Grid.Edges[iE].Mid.z
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    CurlEdge = (CurlUu[iN1] * MetricFV.DualVolume[iN1] + CurlUu[iN2] * MetricFV.DualVolume[iN2]) / 
      (MetricFV.DualVolume[iN1] + MetricFV.DualVolume[iN2])  
    ru[iE] = -TangUu[iE] * (CurlEdge + 2 * Phys.Omega * sinlat) - grad[iE]
    ru[iE] = +TangUu[iE] * (CurlEdge )
  end
  mul!(hE,Inter,Up)
  @. hE *= Uu
  mul!(rpI,Div,hE)
end
function FcnFV4!(F,U,MetricFV,Grid,Cache,Phys)

  pPosS = 1
  pPosEI = Grid.NumFaces
  pPosE = Grid.NumFaces + Grid.NumFacesG
  uPosS = pPosE + 1
  uPosE = pPosE + Grid.NumEdges
  @views rp = F[pPosS:pPosE]
  @views rpI = F[pPosS:pPosEI]
  @views ru = F[uPosS:uPosE]
  @views Up = U[pPosS:pPosE]
  @views UpI = U[pPosS:pPosEI]
  @views Uu = U[uPosS:uPosE]

  CurlUu = Cache.CurlUu
  @views TangUu = Cache.TangUu
  hE = Cache.hE
  grad = Cache.hE
  K = Cache.K
  Curl = Cache.Curl
  Tang = Cache.Tang
  Grad = Cache.Grad
  Div = Cache.Div
  Inter = Cache.Inter

  mul!(TangUu,Tang,Uu)
  KineticEnergy(K,Uu,TangUu,MetricFV,Grid)
  @. K += Up * Phys.Grav
  mul!(CurlUu,Curl,Uu)
  mul!(grad,Grad,K)
  @inbounds for iE = 1 : Grid.NumEdges
    iN1 = Grid.Edges[iE].N[1]
    iN2 = Grid.Edges[iE].N[2]
    x = Grid.Edges[iE].Mid.x
    y = Grid.Edges[iE].Mid.y
    z = Grid.Edges[iE].Mid.z
    sinlat = z / sqrt(x^2 + y^2 + z^2)
    CurlEdge = (CurlUu[iN1] * MetricFV.DualVolume[iN1] + CurlUu[iN2] * MetricFV.DualVolume[iN2]) / 
      (MetricFV.DualVolume[iN1] + MetricFV.DualVolume[iN2])  
    ru[iE] = -TangUu[iE] * (CurlEdge + 2 * Phys.Omega * sinlat) - grad[iE]
    ru[iE] = TangUu[iE] * (2 * Phys.Omega * sinlat) 
  end
  mul!(hE,Inter,Up)
  @. hE *= Uu
  mul!(rpI,Div,hE)
end
