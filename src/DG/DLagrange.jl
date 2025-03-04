function DLagrange(x,xP,iP)
  Df=eltype(x)(0)
  f=eltype(x)(1)
  for i = 1 : size(xP,1)
    if i ≠ iP
      Df = Df * (x-xP[i]) / (xP[iP] - xP[i]) + f / (xP[iP]-xP[i])
      f = f * (x - xP[i]) / (xP[iP] - xP[i])
    end
  end
  return Df
end
