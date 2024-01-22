function Most()

#RiB  = Grav/Tgn*((AbsT-Tgn)*(dz(iz)-z0M))/(VT*VT+Eps)
end
function VertProfile(uf,z0,z,L,tt)
  log(z / z0) - psi(uf, z / L, tt) +
  psi(uf, z0 / L, tt)
end  
function MOSTIteration(uf,z0M,z0H,z,U,theta,thetaS,LandClass,Phys)
  FT = eltype(theta)
  Karm = 0.4
  Pr = 0.74
  Ri_b = min(Phys.Grav / thetaS *(theta - thetaS) * z  / (U * U), 0.2)
  L = (theta - thetaS) / U^2
  for Iter = 1 : 10
    f = Ri_b - Pr * z / L *
      VertProfile(uf, z0H, z, L, HeatTransport()) / 
      VertProfile(uf, z0M, z, L, MomentumTransport())^2 
    LP = L * (FT(1) + sqrt(eps(FT))) + flipsign(sqrt(eps(FT)), L)
    fP = Ri_b - Pr * z / LP *
      VertProfile(uf, z0H, z, LP, HeatTransport()) / 
      VertProfile(uf, z0M, z, LP, MomentumTransport())^2 
    df = (fP - f) /(LP - L)
    LNew = L - f/df
    if abs(LNew-L) < eps(FT)^FT(1/3)
      L  = LNew 
      break
    else
      L = LNew  
    end  
  end  

  CT = Karm^2 / (log(z/z0M) - psi(uf, z / L, MomentumTransport()))^2
  CH = Karm^2 / (log(z/z0M) - psi(uf, z / L, MomentumTransport())) /
    (log(z/z0H) - psi(uf, z / L, HeatTransport()))
  return CT, CH

end

function RichardIteration(z0M,z0H,dz,RiB,LandClass)
  FT = eltype(RiB)
  RiB = max(RiB,-FT(10))
  chi = RiB
  Iter = 0
  Karm = 0.4
  for Iter = 1 : 21
    @show Iter, chi  
    fchi = RiB - RiBulk(dz, z0M, z0H, chi, RiB, LandClass)
    @show fchi
    chiS = chi * (FT(1) + sqrt(eps(FT))) + flipsign(sqrt(eps(FT)), chi)
    fchiS = RiB - RiBulk(dz, z0M, z0H, chiS, RiB, LandClass)
    @show fchiS
    dfdchi = (fchiS - fchi) /(chiS - chi)
    @show fchi, dfdchi
    chiNew = chi - fchi/dfdchi
    @show chi, chiNew
    if abs(chi-chiNew) < eps(FT)^FT(1/3)
      chi  = chiNew 
      break
    elseif Iter>20
      stop
    else
      chi = chiNew  
    end  
  end
  @show Iter
  TermH = log((dz + z0M) / z0H) -      
    PsiH(chi * (FT(1) + z0M / dz), dz, z0M, z0H, RiB, LandClass) +
    PsiH(chi * z0H / dz, dz, z0M, z0H, RiB, LandClass)
  TermM = log((dz+z0M)/z0M) -       
    PsiM(chi*(FT(1) + z0M / dz), dz, z0M, z0H, RiB, LandClass) +
    PsiM(chi * z0M / dz, dz, z0M, z0H, RiB, LandClass)
  Cm = Karm * Karm / TermM / TermM
  Ch = Karm * Karm / TermM / TermH
  return Cm, Ch
end     

