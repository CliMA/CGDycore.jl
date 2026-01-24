mutable struct ModelFEMLin
  RT
  CG
  DG
  pPosS::Int
  pPosE::Int
  uPosS::Int
  uPosE::Int
  Div::AbstractMatrix
  Grad::AbstractMatrix
end

function ModelFEMLin(backend,FTB,RT,CG,DG,Grid,nQuadM,nQuadS::Int,Jacobi)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = pPosE + 1
  uPosE = pPosE + RT.NumG
  DG.M = FEM.MassMatrix(backend,FTB,DG,Grid,nQuadM,Jacobi)
  DG.LUM = lu(DG.M)

  RT.M = FEM.MassMatrix(backend,FTB,RT,Grid,nQuadM,Jacobi)
  RT.LUM = lu(RT.M)
  CG.M = FEM.MassMatrix(backend,FTB,CG,Grid,nQuadM,Jacobi)
  CG.LUM = lu(CG.M)
  Div = FEM.DivMatrix(backend,FTB,RT,DG,Grid,nQuadM,Jacobi)
  Grad = -Div'

  return ModelFEMLin(
    RT,
    CG,
    DG,
    pPosS,
    pPosE,
    uPosS,
    uPosE,
    Div,
    Grad,
  )
end

mutable struct ModelFEMVecI
  RT
  CG
  DG
  pPosS::Int
  pPosE::Int
  uPosS::Int
  uPosE::Int
  Div::AbstractMatrix
  Grad::AbstractMatrix
end


function ModelFEMVecI(backend,FTB,RT,CG,DG,Grid,nQuadM,nQuadS,Jacobi,ExchangeDG,ExchangeCG)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = pPosE + 1
  uPosE = pPosE + RT.NumG
  DG.M = FEM.MassMatrix(backend,FTB,DG,Grid,DG.Order+1,Jacobi)
  DG.LUM = DiagonalMatrix(FTB,DG)
  DE = reshape(DG.LUM.D,1,1,length(DG.LUM.D),1)
  Parallels.ExchangeData3DSendGPU(DE,ExchangeDG)
  Parallels.ExchangeData3DRecvGPU!(DE,ExchangeDG)
# DG.LUM = lu(DG.M)

  RT.M = FEM.MassMatrix(backend,FTB,RT,Grid,nQuadM,Jacobi)
  RT.LUM = BlockDiagonalDualMatrix(FTB,RT,Grid)
# RT.LUM = lu(RT.M)
  Div = FEM.DivMatrix(backend,FTB,RT,DG,Grid,nQuadS,Jacobi)
  Grad = -Div'

  CG.M = FEM.MassMatrix(backend,FTB,CG,Grid,CG.Order,Jacobi)
  CG.LUM = DiagonalMatrix(FTB,CG)
  DE = reshape(CG.LUM.D,1,1,length(CG.LUM.D),1)
  Parallels.ExchangeData3DSendGPU(DE,ExchangeCG)
  Parallels.ExchangeData3DRecvGPU!(DE,ExchangeCG)
# CG.LUM = lu(CG.M)
  return ModelFEMVecI(
    RT,
    CG,
    DG,
    pPosS,
    pPosE,
    uPosS,
    uPosE,
    Div,
    Grad,
  )  
end

mutable struct ModelFEMCons
  ND
  RT
  CG
  DG
  VecDG
  hPosS::Int
  hPosE::Int
  huPosS::Int
  huPosE::Int
  uPosVec::Int
  Div::AbstractMatrix
  Grad::AbstractMatrix
  Curl::AbstractMatrix
  Lapl::AbstractMatrix
end

function ModelFEMCons(backend,FTB,ND,RT,CG,DG,VecDG,Grid,nQuadM,nQuadS,Jacobi)
  hPosS = 1
  hPosE = DG.NumG
  huPosS = hPosE + 1
  huPosE = hPosE + RT.NumG
  uPosVec = VecDG.NumG

  DG.M = FEM.MassMatrix(backend,FTB,DG,Grid,nQuadM,Jacobi)
  DG.LUM = lu(DG.M)

  VecDG.M = FEM.MassMatrix(backend,FTB,VecDG,Grid,nQuadM,Jacobi)
  VecDG.LUM = lu(VecDG.M)

  RT.M = FEM.MassMatrix(backend,FTB,RT,Grid,nQuadM,Jacobi)
  RT.LUM = lu(RT.M)
  Div = FEM.DivMatrix(backend,FTB,RT,DG,Grid,nQuadS,Jacobi)
  Grad = -Div'

  ND.M = FEM.MassMatrix(backend,FTB,ND,Grid,nQuadM,Jacobi)
  ND.LUM = lu(ND.M)
  CG.M = FEM.MassMatrix(backend,FTB,CG,Grid,nQuadM,Jacobi)
  CG.LUM = lu(CG.M)
  Curl = FEM.CurlMatrix(backend,FTB,ND,DG,Grid,nQuadS,Jacobi)
  Lapl = spzeros(0,0)
  return ModelFEMCons(
    ND,
    RT,
    CG,
    DG,
    VecDG,
    hPosS,
    hPosE,
    huPosS,
    huPosE,
    uPosVec,
    Div,
    Grad,
    Curl,
    Lapl,
  )  
end
