mutable struct ModelFEM
  ND::HCurlElement
  RT::HDivElement
  DG::ScalarElement
  pPosS::Int
  pPosE::Int
  uPosS::Int
  uPosE::Int
end

function ModelFEM(ND,RT,DG)
  pPosS = 1
  pPosE = DG.NumG
  uPosS = pPosE + 1
  uPosE = pPosE + RT.NumG
  return ModelFEM(
    ND,
    RT,
    DG,
    pPosS,
    pPosE,
    uPosS,
    uPosE,
  )  
end
