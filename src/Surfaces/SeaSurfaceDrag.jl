function DragCoeffSea(U,zPL,TSurf,Phys)
  
    RhoLoc = U[1]
    RhoThLoc = U[5]
    ThLoc = RhoThLoc / RhoLoc
    VT = U[2]^2 + U[3]^2
    One = 1.0
    Two = 2.0
#   Turbulent Prandtl Number (after Businger)
    R = 7.4e-1
  
    
    Fpot      = (Phys.Rd*RhoThLoc/Phys.p0)^(-Phys.kappa/(One-Phys.kappa))
    ThPotSurf = TSurf*Fpot
  
#   First Estimate for Bulk Richardson Number RiB0 = RiBT
    RiBT      = Phys.Grav*zPL*(ThLoc-ThPotSurf)/(ThLoc*(VT^Two))

    if RiBT >0.25   # critical Richardson number: 0.25
#   Stable Conditions
      FM = One/((One+4.7d0*RiBT)^2.0d0)
      FH = FM
    else 
#   Unstable Conditions
      cM = 7.4d0*(Karm^2.0d0)*9.4d0*sqrt(zPL/(zRauh1))/((log(zPL/zRauh1))^2.0d0)
      cH = 5.3d0*cM/7.4d0
      FM = One-9.4d0*RiBT/(One+cM*sqrt(ABS(RiBT)))
      FH = One-9.4d0*RiBT/(One+cH*sqrt(ABS(RiBT)))
    end

    for it = 0 : 5      

#     Friction velocity for momentum
      ustar  = sqrt(FM*(Karm*VT/(log(zPL/zRauh1)))^Two)

      if RiBT > 0.25
        wstar=Zero
      else
        wstar = sqrt((Karm^Two/(R*log(zPl/zRauh1)*log(zpl/zRauhT))*VT^3*ABS(RiBT)*FH)^(3/2))
        ustar = max(ustar,wstar)
      end

!     New roughness lenght for momentum (after Charnock 1955,Nikuradse 1933)
      zRauh2=1.1d-2*ustar^Two/Grav
      zRauh2= max(zRauh2,1d-5)

      zRauhT=zRauh1
      zRauhT=min(zRauhT,0.1d0)

      if (ABS(zRauh1-zRauh2)<=10e-6*ABS(zRauh2).AND.zRauh1>1d-5)   
        zRauh1=zRauh2
        exit
      eeSE
        zRauh1=zRauh2
      end
    end

!   Friction velocity for momentum        
    ustar  = sqrt(FM*(Karm*VT/(log(zPL/zRauh1)))^Two)
    U10  = log(10d0/zRauh1)*ustar/(sqrt(FM)*Karm)


!   Drag coefficient for momentum, heat and moisture
    DragM  = (ustar/(VT))^Two
    DragH  = (Karm/(log(zPL/zRauh1)))^Two*FH/(R*(log(zRauh1/zRauhT)*FH/(log(zPL/zRauh1)*sqrt(FM))+One))          

end
