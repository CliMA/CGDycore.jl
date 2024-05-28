function Fcn!(F,U,Model)

  @views pU = U[Model.pPosS:Model.pPosE]
  @views uU = U[Model.uPosS:Model.uPosE]
  @views pF = F[Model.pPosS:Model.pPosE]
  @views uF = F[Model.uPosS:Model.uPosE]


end
