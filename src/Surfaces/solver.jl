function Most()

#RiB  = Grav/Tgn*((AbsT-Tgn)*(dz(iz)-zRauh))/(VT*VT+Eps)
end
function RichardIteration(zRauh,zRauhT,dz,RiB)
  FT = eltype(RiB)
  RiB = MAX(RiB,-FT(10))
  chi = RiB
  Iter=0
  for Iter = 1 : 21
    fchi = FixPoint(chi,RiB)
    chiS = chi*(FT(1)+sqrt(eps(FT)))+sign(sqrt(eps(FT)),chi)
    fchiS = FixPoint(chiSi,RiB)
    dfdchi = (fchiS - fchi) /(chiS - chi)
    chiNew = chi - fchi/dfdchi
    if abs(chi-chiNew) < eps(FT)^FT(1/3)
      break
    elseif Iter>20
      stop
    else
      chi = chiNew  
    end  
  end
  LandClass = 1
  TermH = log((dz+zRauh)/zRauhT) -      
    PsiH(chiNew*(FT(1)One+zRauh/dz),dz,zRauh,zRauhT,RiB,LandClass) +
    PsiH(chiNew*zRauhT/dz,dz,zRauh,zRauhT,RiB,LandCLass)
  TermM = log((dz+zRauh)/zRauh) -       
    PsiM(chiNew*(FT(1) + zRauh/dz),dz,zRauh,zRauhT,RiB,LandCLass) +
    PsiM(chiNew*zRauh/dz,dz,zRauh,zRauhT,RiB,LandCLass)
  Cm = Karm * Karm / TermM / TermM
  Ch = Karm Karm / TermM / TermH
end     

@inline function FixPoint(chi,RiB)
  RiB - RiBulk(chi,RiB)
end  
