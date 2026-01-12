function FcnLinShallow!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  @. F = 0
  DG = Model.DG
  RT = Model.RT

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]

  DivRhs!(backend,FTB,Fh,DG,Uu,RT,Grid,Grids.Quad(),QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fh)

  CrossRhs!(backend,FTB,Fu,RT,Uu,RT,Grid,Grids.Quad(),QuadOrdS,Jacobi!)
  GradRhs!(backend,FTB,Fu,RT,Uh,DG,Grid,Grids.Quad(),QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)

end

function FcnHeat!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  DG = Model.DG
  Lapl = Model.Lapl

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]

  mul!(Fh,Lapl,Uh)
  ldiv!(DG.LUM,Fh)
end

function FcnVecINonLinShallow!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  @. F = 0
  DG = Model.DG
  CG = Model.CG
  RT = Model.RT
  Grav =  9.80616

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views UCachep = UCache[Model.pPosS:Model.pPosE]
  @views k = UCache[Model.pPosS:Model.pPosE]
  @views q = UCache[Model.uPosE+1:end]
  @views hUu = UCache[Model.uPosS:Model.uPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]

  CurlVel!(q,CG,Uu,RT,QuadOrdS,Grid.Type,Grid,Jacobi)
  CrossRhs!(backend,FTB,Fu,RT,q,CG,Uu,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  KineticEnergy!(backend,FTB,k,DG,Uu,RT,Grid,Grids.Quad(),QuadOrdM,Jacobi)
  @. k += Grav * Uh
  GradRhs!(backend,FTB,Fu,RT,k,DG,Grid,Grids.Quad(),QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)
  ProjecthScalaruHDivHDiv!(backend,FTB,hUu,RT,Uh,DG,Uu,RT,Grid,RT.Type,QuadOrdM,Jacobi)
  DivRhs!(backend,FTB,Fh,DG,hUu,RT,Grid,DG.Type,QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fh)

#=
  views Uh = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views UCachep = UCache[Model.pPosS:Model.pPosE]
  @views k = UCache[Model.pPosS:Model.pPosE]
  @views UCacheu = UCache[Model.uPosS:Model.uPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]
  @views UCacheq = UCache[Model.uPosE+1:end]

  CurlVel!(UCacheq,CG,Uu,RT,QuadOrdS,Grid.Type,Grid,Jacobi)
  CrossRhs!(backend,FTB,Fu,UCacheq,CG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  GradKinHeight!(backend,FTB,Fu,Uh,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)
  ProjecthScalaruHDivHDiv!(backend,FTB,UCacheu,RT,Uh,DG,Uu,RT,Grid,RT.Type,QuadOrdM,Jacobi)
  DivRhs!(backend,FTB,Fh,DG,UCacheu,RT,Grid,DG.Type,QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fh)
=#  
end

function FcnVecINonLinShallowRT!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

  @. F = 0
  DG = Model.DG
  CG = Model.CG
  RT = Model.RT
  ND = Model.ND
  Div = Model.Div
  Grad = Model.Grad

  @views Uh = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views UCachep = UCache[Model.pPosS:Model.pPosE]
  @views k = UCache[Model.pPosS:Model.pPosE]
  @views UCacheu = UCache[Model.uPosS:Model.uPosE]
  @views Fh = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]
  @views UCacheq = UCache[Model.uPosE+1:end]

  CurlVel!(UCacheq,CG,Uu,RT,QuadOrdS,Grid.Type,Grid,Jacobi)
  CrossRhs!(backend,FTB,Fu,UCacheq,CG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  GradKinHeight!(backend,FTB,Fu,Uh,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  ldiv!(RT.LUM,Fu)
  ProjecthScalaruHDivHDiv!(backend,FTB,UCacheu,RT,Uh,DG,Uu,RT,Grid,RT.Type,QuadOrdM,Jacobi)
  DivRhs!(backend,FTB,Fh,DG,UCacheu,RT,Grid,DG.Type,QuadOrdS,Jacobi)
  ldiv!(DG.LUM,Fh)
end

function FcnConsNonLinShallow!(backend,FTB,F,U,Model,Grid,nQuad,Jacobi;UCache=nothing)
  @. F = 0
  hPosS = Model.hPosS
  hPosE = Model.hPosE
  huPosS = Model.huPosS
  huPosE = Model.huPosE
  uPosVec = Model.uPosVec

  @views Uh = U[hPosS:hPosE]
  @views Uhu = U[huPosS:huPosE]
  @views Fh = F[hPosS:hPosE]
  @views Fhu = F[huPosS:huPosE]
  uRec = zeros(FTB, uPosVec)

  # Tendency h
  DivRhs!(backend,FTB,Fh,Model.DG,Uhu,Model.RT,Grid,Grid.Type,nQuad,Jacobi)
  ldiv!(Model.DG.LUM,Fh)
  # Tendency hu
  InterpolateScalarHDivVecDG!(backend,FTB,uRec,Model.VecDG,Uh,Model.DG,Uhu,Model.RT,Grid,Grid.Type,nQuad,Jacobi)
  DivMomentumVector!(backend,FTB,Fhu,Model.RT,Uhu,Model.RT,uRec,Model.VecDG,Grid,Grid.Type,nQuad,Jacobi)
  if Grid.Form == "Sphere"
    CrossRhs!(backend,FTB,Fhu,Model.RT,Uhu,Model.RT,Grid,Grid.Type,nQuad,Jacobi)
  end
  GradHeightSquared!(backend,FTB,Fhu,Model.RT,Uh,Model.DG,Grid,Grid.Type,nQuad,Jacobi)
  ldiv!(Model.RT.LUM,Fhu)
end

function FcnNonLinShallowKent!(backend,FTB,F,U,Model,Grid,QuadOrdM,QuadOrdS,Jacobi;UCache)

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
  @views UCacheq = UCache[Model.uPosE+1:end]

  ProjectHDivHCurl!(backend,FTB,UCacheu,ND,Uu,RT,
    Grid,RT.Type,QuadOrdM,Jacobi)
  mul!(UCachep,Curl,UCacheu)
  ldiv!(DG.LUM,UCachep)
  CrossRhs!(backend,FTB,Fu,UCacheq,CG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
  GradKinHeight!(backend,FTB,Fu,Uh,DG,Uu,RT,RT,Grid,RT.Type,QuadOrdS,Jacobi)
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
