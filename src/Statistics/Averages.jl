function AverageInTime!(UAver,U,Iter)

  @. UAver = UAver + (U - UAver) / Iter
end
