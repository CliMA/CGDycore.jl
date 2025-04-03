function FcnLinShallow!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  DG = Model.DG
  RT = Model.RT
  Div = Model.Div
  Grad = Model.Grad

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]
  @views Gravh = UCache[Model.pPosS:Model.pPosE]

  @. Gravh = 9.80616 * Uh
  GradRhs!(backend,FTB,Fu,Gravh,DG,RT,Grid,Grid.Type,QuadOrdS,Jacobi)
  CrossRhs!(backend,FTB,Fu,RT,Uu,RT,Grid,Grid.Type,QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)

  DivRhs!(backend,FTB,Fh,DG,Uu,RT,Grid,Grid.Type,QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fh)
end

function FcnHeat!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  DG = Model.DG
  Lapl = Model.Lapl

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]

  mul!(Fh,Lapl,Uh)
  ldiv!(DG.LUM,Fh)
end

function FcnNonLinShallow!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  @. F = 0
  DG = Model.DG
  RT = Model.RT
  ND = Model.ND
  Div = Model.Div
  Grad = Model.Grad
  Curl = Model.Curl

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views UCachep = UCache[Model.pPosS:Model.pPosE]
  @views k = UCache[Model.pPosS:Model.pPosE]
  @views UCacheu = UCache[Model.uPosS:Model.uPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]

  ProjectHDivHCurl!(backend,FTB,UCacheu,ND,Uu,RT,
    Grid,RT.Type,QuadOrdM,Jacobi)
  mul!(UCachep,Curl,UCacheu)
  ldiv!(DG.LUM,UCachep)
# CurlVel!(UCachep,DG,Uu,RT,QuadOrdS,Grid.Type,Grid,Jacobi)
  CrossRhs!(backend,FTB,Fu,UCachep,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  GradKinHeight!(backend,FTB,Fu,Uh,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
# ProjectKE!(backend,FTB,k,DG,Uu,RT,Grid.Type,Grid,QuadOrdS,Jacobi)
# @. k += 9.81 * Uh
# GradRhs!(backend,FTB,Fu,k,DG,RT,Grid,Grid.Type,QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)
  ProjecthScalaruHDivHDiv!(backend,FTB,UCacheu,RT,Uh,DG,Uu,RT,Grid,RT.Type,QuadOrdM,Jacobi)
  DivRhs!(backend,FTB,Fh,DG,UCacheu,RT,Grid,DG.Type,QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fh)
end

function Curl!(backend,FTB,uCurl,DG,Uu,RT,ND,Grid,Jacobi,QuadOrdM,Curl,UCacheu)
  ProjectHDivHCurl!(backend,FTB,UCacheu,ND,Uu,RT,
    Grid,RT.Type,QuadOrdM,Jacobi)
  mul!(uCurl,Curl,UCacheu)
  ldiv!(DG.LUM,uCurl)
end
