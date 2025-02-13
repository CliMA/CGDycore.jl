function FcnLinShallow!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  DG = Model.DG
  RT = Model.RT
  Div = Model.Div
  Grad = Model.Grad

  @views Up = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views Fp = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]

  mul!(Fu,Grad,Up)
  ldiv!(RT.LUM,Fu)

  mul!(Fp,Div,Uu)
  ldiv!(DG.LUM,Fp)
end

function FcnHeat!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  DG = Model.DG
  Lapl = Model.Lapl

  @views Up = U[Model.pPosS:Model.pPosE]
  @views Fp = F[Model.pPosS:Model.pPosE]

  mul!(Fp,Lapl,Up)
  ldiv!(DG.LUM,Fp)
end

function FcnNonLinShallow!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  @. F = 0
  DG = Model.DG
  RT = Model.RT
  ND = Model.ND
  Div = Model.Div
  Grad = Model.Grad
  Curl = Model.Curl

  @views Up = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views UCachep = UCache[Model.pPosS:Model.pPosE]
  @views k = UCache[Model.pPosS:Model.pPosE]
  @views UCacheu = UCache[Model.uPosS:Model.uPosE]
  @views Fp = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]

  ProjectHDivHCurl!(backend,FTB,UCacheu,ND,Uu,RT,
    Grid,RT.Type,QuadOrdM,Jacobi)
  mul!(UCachep,Curl,UCacheu)
  ldiv!(DG.LUM,UCachep)


# CurlVel(UCachep,DG,Uu,RT,QuadOrdS,Grid.Type,Grid,Jacobi)


  CrossRhs!(backend,FTB,Fu,UCachep,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  GradKinHeight!(backend,FTB,Fu,Up,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)

  ProjecthScalaruHDivHDiv!(backend,FTB,UCacheu,RT,Up,DG,Uu,RT,Grid,RT.Type,QuadOrdM,Jacobi)
  DivRhs!(backend,FTB,Fp,DG,UCacheu,RT,Grid,DG.Type,QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fp)

end

function Curl!(backend,FTB,uCurl,DG,Uu,RT,ND,Grid,Jacobi,QuadOrdM,Curl,UCacheu)
  ProjectHDivHCurl!(backend,FTB,UCacheu,ND,Uu,RT,
    Grid,RT.Type,QuadOrdM,Jacobi)
  mul!(uCurl,Curl,UCacheu)
  ldiv!(DG.LUM,uCurl)
end
