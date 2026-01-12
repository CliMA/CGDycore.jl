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

function ModelFEMLin(backend,FTB,RT,CG,DG,Grid,nQuadM,nQuadS,Jacobi)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = pPosE + 1
  uPosE = pPosE + RT.NumG
  @show nQuadM
  DG.M = FEM.MassMatrix(backend,FTB,DG,Grid,nQuadM,Jacobi)
  DG.LUM = lu(DG.M)

  RT.M = FEM.MassMatrix(backend,FTB,RT,Grid,nQuadM,Jacobi)
  RT.LUM = lu(RT.M)
  CG.M = FEM.MassMatrix(backend,FTB,CG,Grid,nQuadM,Jacobi)
  CG.LUM = lu(CG.M)
  Div = FEM.DivMatrix(backend,FTB,RT,DG,Grid,nQuadS,Jacobi)
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
  ND
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


function ModelFEMVecI(backend,FTB,ND,RT,CG,DG,Grid,nQuadM,nQuadS,Jacobi)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = pPosE + 1
  uPosE = pPosE + RT.NumG
  DG.M = FEM.MassMatrix(backend,FTB,DG,Grid,nQuadM,Jacobi)
  DG.LUM = lu(DG.M)

  RT.M = FEM.MassMatrix(backend,FTB,RT,Grid,nQuadM,Jacobi)
  RT.LUM = lu(RT.M)
  Div = FEM.DivMatrix(backend,FTB,RT,DG,Grid,nQuadS,Jacobi)
  Grad = -Div'

  ND.M = FEM.MassMatrix(backend,FTB,ND,Grid,nQuadM,Jacobi)
  ND.LUM = lu(ND.M)
  CG.M = FEM.MassMatrix(backend,FTB,CG,Grid,nQuadM,Jacobi)
  CG.LUM = lu(CG.M)
  return ModelFEMVecI(
    ND,
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
