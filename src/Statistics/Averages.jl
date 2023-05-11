function AverageInTime!(UAver,U,Iter)

  @. UAver = UAver + (U - UAver) / Iter
  Iter += 1
end
