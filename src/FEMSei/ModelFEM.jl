mutable struct ModelFEM
  ND
  RT
  DG
  pPosS::Int
  pPosE::Int
  uPosS::Int
  uPosE::Int
  Div::AbstractMatrix
  Grad::AbstractMatrix
  Curl::AbstractMatrix
  Lapl::AbstractMatrix
end

function ModelFEM(backend,FTB,DG,Grid,nQuad,Jacobi)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = 0
  uPosE = 0
  ND = nothing
  RT = nothing
  Div = spzeros(0,0)
  Grad = spzeros(0,0)
  Curl = spzeros(0,0)
  Lapl = FEMSei.LaplMatrix(backend,FTB,DG,DG,Grid,nQuad,Jacobi)
  DG.M = FEMSei.MassMatrix(backend,FTB,DG,Grid,nQuad,Jacobi)
  DG.LUM = lu(DG.M)
  return ModelFEM(
    ND,
    RT,
    DG,
    pPosS,
    pPosE,
    uPosS,
    uPosE,
    Div,
    Grad,
    Curl,
    Lapl,
  )  
end

function ModelFEM(backend,FTB,ND,RT,DG,Grid,nQuad,Jacobi)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = pPosE + 1
  uPosE = pPosE + RT.NumG
  DG.M = FEMSei.MassMatrix(backend,FTB,DG,Grid,nQuad,Jacobi)
  DG.LUM = lu(DG.M)

  RT.M = FEMSei.MassMatrix(backend,FTB,RT,Grid,nQuad,Jacobi)
  RT.LUM = lu(RT.M)
  Div = FEMSei.DivMatrix(backend,FTB,RT,DG,Grid,nQuad,Jacobi)
  Grad = -Div'

  ND.M = FEMSei.MassMatrix(backend,FTB,ND,Grid,nQuad,Jacobi)
  ND.LUM = lu(ND.M)
  Curl = FEMSei.CurlMatrix(backend,FTB,ND,DG,Grid,nQuad,Jacobi)
  Lapl = spzeros(0,0)
  return ModelFEM(
    ND,
    RT,
    DG,
    pPosS,
    pPosE,
    uPosS,
    uPosE,
    Div,
    Grad,
    Curl,
    Lapl,
  )  
end
