function RiBulk(z,z0M,z0h,Chi,RiB,LandClass)
  FT = eltype(z)
  PsiHeat = log((z+z0)/z0h) -
    PsiH(Chi*(FT(1)+z0/z),z,z0M,z0h,RiB,LandClass) +
    PsiH(Chi*z0h/z,z,z0M,z0h,RiB,LandClass)
  PsiMom = log((z+z0M)/z0) -
    PsiM(Chi*(FT(1)+z0M/z),z,z0M,z0h,RiB,LandClass) +
    PsiM(Chi*z0M/z,z,z0M,z0h,RiB,LandClass)
  return Chi * PsiHeat / (PsiMom * PsiMom)
end

function PsiH(Chi,z,z0M,z0h,RiB,LandClass)
  if LandClass <= 8
    a = FT(-5.3)
    b = FT(1.1)
    if RiB >= FT(0)
      PsiH  = a*log(Chi+(FT(1)+ Chi^b)^(FT(1)/b))
    else
      PsiH  = (PsiHK(Chi)+Chi*Chi*PsiH_C(Chi))/(FT(1)+Chi*Chi)
    end
    PsiH = min(FT(0.9)*log((z+z0h)/z0M),PsiH)
  elseif LandClass == 10
    PsiH = PsiH_SHEBA(Chi)
  end
  return PsiH
end

function PsiM(Chi,z,z0M,z0h,RiB,LandClass)
  FT = eltype(z)
  if LandClass <= 8
    a = FT(-6.1)
    b =  FT(2.5)
    if RiB >= FT(0) 
      PsiM  = a*log(Chi+(FT(1)+Chi**b)**(One/b))
    else
      PsiM  = (PsiMK(Chi)+Chi*Chi*PsiM_C(Chi))/(FT(1)+Chi*Chi)
    end
    PsiM = min(FT(0.9)*log((z+z0h)/z0M),PsiM)
  elseif LandClass == 10
    PsiM = PsiM_SHEBA(Chi)
  end   
  return PsiM
end   

function PsiHK(Chi)
  FT = eltype(Chi)
  alpha = FT(16)
  Var1  = (FT(1)-alpha*Chi)^FT(1/4)
  return FT(2) * log((FT(1) + Var1 * Var1) / FT(2))
end

function PsiMK(Chi)
  FT = eltype(Chi)
  alpha = FT(16)
  Var1  = (FT(1)-alpha*Chi)^FT(1/4)
  return FT(2) * log((FT(1) + Var1) / FT(2))  +
    log((FT(1) + Var1 * Var1) / FT(2)) -
    FT(2) * atan(Var1) + FT(0.5) * pi
end

function PsiH_C(Chi)
  FT = eltype(Chi)
  beta = 34.d0
  Var2 = (FT(1)-beta*Chi)**(FT(1)/FT(3))
  PsiH_C = FT(3)/FT(2)*log((Var2**FT(2)+Var2+FT(1))/FT(3)) -
    sqrt(FT(3))*atan((FT(2)*Var2+FT(1))/sqrt(FT(3))) +
    Pi/sqrt(FT(3))
end

function PsiM_C(Chi)
  beta = 10.d0
  Var2 = (FT(1)-beta*Chi)**(FT(1)/FT(3))
  PsiM_C = FT(3)/FT(2)*log((Var2**FT(2)+Var2+FT(1))/FT(3)) -
    sqrt(FT(3))*atan((FT(2)*Var2+FT(1))/sqrt(FT(3))) +
    pi/sqrt(FT(3))
end

# formula from SHEBA flux–profile relationships in the stable
# atmospheric boundary layer by Gravech et al. 2007
function PsiH_SHEBA(Chi)
  FT = eltype(Chi)

  ah = 5.d0
  bh = 5.d0
  ch = 3.d0
  grbh = sqrt(ch^FT(2)-FT(4))

  PsiH_SHEBA = - bh/FT(2) * log(FT(1)+ch*Chi+Chi) +
    (-ah/grbh + bh*ch/(FT(2)*grbh)) *
    (log((FT(2)*Chi+ch-grbh)/(FT(2)*Chi+ch+grbh)) -
    log((ch-grbh)/(ch+grbh)))

end

function PsiM_SHEBA(Chi)
  FT = eltype(Chi)

  am = 5.d0
  bm = am/6.5d0
  grbm = ((FT(1)-bm)/bm)^(FT(1)/FT(3))
  Var3 = (FT(1)+Chi)^(FT(1)/FT(3))

  PsiM_SHEBA = -FT(3)*am/bm*(Var3-FT(1))+(am*grbm)/(FT(2)*bm) *
    (FT(2)*log((Var3+grbm)/(FT(1)+grbm)) -
    log(Var3^FT(2)-Var3*grbm+grbm^FT(2))/(FT(1)-grbm+grbm^FT(2)) +
    FT(2)*sqrt(FT(3))*(atan((FT(2)*Var3-grbm)/(sqrt(FT(3))*grbm)) -
    atan((FT(2)-grbm)/(sqrt(FT(3))*grbm))))
end



