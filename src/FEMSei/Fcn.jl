function Fcn!(F,U,Model)

  DG = Model.DG
  RT = Model.RT
  Div = Model.Div
  Grad = Model.Grad

  @views Up = U[Model.pPosS:Model.pPosE]
  @views Uu = U[Model.uPosS:Model.uPosE]
  @views Fp = F[Model.pPosS:Model.pPosE]
  @views Fu = F[Model.uPosS:Model.uPosE]

  mul!(Fp,Div,Uu)
  ldiv!(DG.LUM,Fp)
  mul!(Fu,Grad,Up)
  ldiv!(RT.LUM,Fu)
end
