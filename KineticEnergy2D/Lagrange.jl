function Lagrange(x,xP,iP)
  f=1
  for i=1:size(xP,1)
     if i!=iP
        f=f*(x-xP[i])/(xP[iP]-xP[i])
     end
  end
  return f
end


function DLagrange(x,xP,iP)
  Df=0.0
  f=1.0
  for i=1:size(xP,1)
    if i!=iP
      Df=Df*(x-xP[i])/(xP[iP]-xP[i])+f/(xP[iP]-xP[i])
      f=f*(x-xP[i])/(xP[iP]-xP[i])
    end
  end
  return Df
end
