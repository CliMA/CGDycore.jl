function GlobalIntegral(c,CG,Global)
  SumLoc = sum(c.*CG.MMass)
  SumGlobal = MPI.Allreduce(SumLoc, +, MPI.COMM_WORLD)
end
