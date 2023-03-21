function CoriolisColumn!(FuC,FvC,uC,vC,RhoC,Fe,X,J,Omega)
  Nz = size(FuC,3)
  OrdPoly = Fe.OrdPoly

  @inbounds for iz = 1 : Nz 
    @inbounds for i = 1 : OrdPoly + 1
      @inbounds for j = 1 : OrdPoly + 1
        x = 0.5 * (X[i,j,1,1,iz] + X[i,j,2,1,iz])  
        y = 0.5 * (X[i,j,1,2,iz] + X[i,j,2,2,iz])  
        z = 0.5 * (X[i,j,1,3,iz] + X[i,j,2,3,iz])  
#       lon,lat,r = cart2sphere(x,y,z)
        r = sqrt(x^2 + y^2 + z^2);
        sinlat = z/r
        W = -2.0 * Omega * sinlat * (J[i,j,1,iz] + J[i,j,2,iz])
        FuC[i,j,iz] -= RhoC[i,j,iz] * vC[i,j,iz] * W
        FvC[i,j,iz] += RhoC[i,j,iz] * uC[i,j,iz] * W
      end 
    end 
  end 
end

function CurlColumn!(FuC,FvC,Fw,uC,vC,w,RhoC,Fe,dXdxI,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4, TCacheCol1, TCacheCol2, 
    TCacheCol3, TCacheCCC1, TCacheCCC2, TCacheCCC3, TCacheCCC4 = ThreadCache
  Nz = size(FuC,3)
  OrdPoly = Fe.OrdPoly
  D = Fe.DS

  @views tempuZ = TCacheCol1[Threads.threadid()]
  @views tempvZ = TCacheCol2[Threads.threadid()]
  @views tempwZ = TCacheCol3[Threads.threadid()]
  @views U = TCacheCCC1[Threads.threadid()]
  @views V = TCacheCCC2[Threads.threadid()]
  @views W = TCacheCCC3[Threads.threadid()]
  @views temp = TCacheCCC4[Threads.threadid()]

  @views FluxUZ = TCacheC2[Threads.threadid()]
  @views FluxVZ = TCacheC3[Threads.threadid()]
  @views FluxWZ = TCacheC4[Threads.threadid()]

  @inbounds for iz = 1 : Nz 
    # Dz*(dx33*v - dx32*w)
    @views @. tempuZ[:,:,1,iz] = dXdxI[:,:,1,iz,3,3] * vC[:,:,iz] - dXdxI[:,:,1,iz,3,2] * w[:,:,iz]
    @views @. tempuZ[:,:,2,iz] = dXdxI[:,:,2,iz,3,3] * vC[:,:,iz] - dXdxI[:,:,2,iz,3,2] * w[:,:,iz+1]
    # Dz*(dx33*u - dx31*w)
    @views @. tempvZ[:,:,1,iz] = dXdxI[:,:,1,iz,3,3] * uC[:,:,iz] - dXdxI[:,:,1,iz,3,1] * w[:,:,iz]
    @views @. tempvZ[:,:,2,iz] = dXdxI[:,:,2,iz,3,3] * uC[:,:,iz] - dXdxI[:,:,2,iz,3,1] * w[:,:,iz+1]
    # Dz*(dx32*u - dx31*v)
    @views @. tempwZ[:,:,1,iz] = dXdxI[:,:,1,iz,3,2] * uC[:,:,iz] - dXdxI[:,:,1,iz,3,1] * vC[:,:,iz]
    @views @. tempwZ[:,:,2,iz] = dXdxI[:,:,2,iz,3,2] * uC[:,:,iz] - dXdxI[:,:,2,iz,3,1] * vC[:,:,iz] 
  end  
  @inbounds for iz = 2 : Nz 
    @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz] - tempuZ[:,:,2,iz-1])
    @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz] - tempvZ[:,:,2,iz-1])
    @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz] - tempwZ[:,:,2,iz-1])
    @views @. FuC[:,:,iz] += RhoC[:,:,iz] * (-vC[:,:,iz] * FluxWZ - w[:,:,iz] * FluxVZ)
    @views @. FvC[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz] * FluxWZ - w[:,:,iz] * FluxUZ)
    @views @. Fw[:,:,iz] += RhoC[:,:,iz] * ( uC[:,:,iz] * FluxVZ + vC[:,:,iz] * FluxUZ)
  end
  @inbounds for iz = 1 : Nz - 1 
    @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz+1] - tempuZ[:,:,2,iz])
    @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz+1] - tempvZ[:,:,2,iz])
    @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz+1] - tempwZ[:,:,2,iz])
    @views @. FuC[:,:,iz] += RhoC[:,:,iz] * (-vC[:,:,iz] * FluxWZ - w[:,:,iz+1] * FluxVZ)
    @views @. FvC[:,:,iz] += RhoC[:,:,iz] * ( uC[:,:,iz] * FluxWZ - w[:,:,iz+1] * FluxUZ)
    @views @. Fw[:,:,iz+1] += RhoC[:,:,iz] * (uC[:,:,iz] * FluxVZ + vC[:,:,iz] * FluxUZ)
  end    
  @inbounds for iz = 1 : Nz 
    @. U = 0.0
    @. V = 0.0
    @. W = 0.0

# uDot = - v*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))  
#        - w*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w))
# vDot =   u*(Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v))
#        - w*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
# wDot =   u*(Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)) 
#          v*(Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w))
#   U = Dx*(dx13*v - dx12*w) + Dy*(dx23*v - dx22*w) + Dz*(dx33*v - dx32*w)    
#   V = Dx*(dx13*u - dx11*w) + Dy*(dx23*u - dx21*w) + Dz*(dx33*u - dx31*w)
#   W = Dx*(dx12*u - dx11*v) + Dy*(dx22*u - dx21*v) + Dz*(dx32*u - dx31*v)
    @views @. temp[:,:,1] = dXdxI[:,:,1,iz,1,3] * vC[:,:,iz] - dXdxI[:,:,1,iz,1,2] * w[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,2,iz,1,3] * vC[:,:,iz] - dXdxI[:,:,2,iz,1,2] * w[:,:,iz+1]
    DerivativeX!(U,temp,D)
    @views @. temp[:,:,1] = dXdxI[:,:,1,iz,2,3] * vC[:,:,iz] - dXdxI[:,:,1,iz,2,2] * w[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,2,iz,2,3] * vC[:,:,iz] - dXdxI[:,:,2,iz,2,2] * w[:,:,iz+1]
    DerivativeY!(U,temp,D)
    @views @. temp[:,:,1] = dXdxI[:,:,1,iz,1,3] * uC[:,:,iz] - dXdxI[:,:,1,iz,1,1] * w[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,2,iz,1,3] * uC[:,:,iz] - dXdxI[:,:,2,iz,1,1] * w[:,:,iz+1]
    DerivativeX!(V,temp,D)
    @views @. temp[:,:,1] = dXdxI[:,:,1,iz,2,3] * uC[:,:,iz] - dXdxI[:,:,1,iz,2,1] * w[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,2,iz,2,3] * uC[:,:,iz] - dXdxI[:,:,2,iz,2,1] * w[:,:,iz+1]
    DerivativeY!(V,temp,D)
    @views @. temp[:,:,1] = dXdxI[:,:,1,iz,1,2] * uC[:,:,iz] - dXdxI[:,:,1,iz,1,1] * vC[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,2,iz,1,2] * uC[:,:,iz] - dXdxI[:,:,2,iz,1,1] * vC[:,:,iz]
    DerivativeX!(W,temp,D)
    @views @. temp[:,:,1] = dXdxI[:,:,1,iz,2,2] * uC[:,:,iz] - dXdxI[:,:,1,iz,2,1] * vC[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,2,iz,2,2] * uC[:,:,iz] - dXdxI[:,:,2,iz,2,1] * vC[:,:,iz]
    DerivativeY!(W,temp,D)

    @views @. U[:,:,1] += 0.5 * (tempuZ[:,:,2,iz] - tempuZ[:,:,1,iz])
    @views @. U[:,:,2] += 0.5 * (tempuZ[:,:,2,iz] - tempuZ[:,:,1,iz])
    @views @. V[:,:,1] += 0.5 * (tempvZ[:,:,2,iz] - tempvZ[:,:,1,iz])
    @views @. V[:,:,2] += 0.5 * (tempvZ[:,:,2,iz] - tempvZ[:,:,1,iz])
    @views @. W[:,:,1] += 0.5 * (tempwZ[:,:,2,iz] - tempwZ[:,:,1,iz])
    @views @. W[:,:,2] += 0.5 * (tempwZ[:,:,2,iz] - tempwZ[:,:,1,iz])

    @views @. FuC[:,:,iz] += RhoC[:,:,iz] * (-vC[:,:,iz] * W[:,:,1] - w[:,:,iz] * V[:,:,1] -
      vC[:,:,iz] * W[:,:,2] - w[:,:,iz+1] * V[:,:,2]) 
    @views @. FvC[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz] * W[:,:,1] - w[:,:,iz] * U[:,:,1] +
      uC[:,:,iz] * W[:,:,2] - w[:,:,iz+1] * U[:,:,2]) 
    @views @. Fw[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz] * V[:,:,1] + vC[:,:,iz] * U[:,:,1])
    @views @. Fw[:,:,iz+1] += RhoC[:,:,iz] * (uC[:,:,iz] * V[:,:,2] + vC[:,:,iz] * U[:,:,2])
  end 
end

function RhoGradColumn!(FuC,FvC,Fw,pC,RhoC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2 = ThreadCache
  Nz = size(FuC,3)
  OrdPoly = Fe.OrdPoly
  D = Fe.DS

  @views DXpC = TCacheC1[Threads.threadid()]
  @views DYpC = TCacheC2[Threads.threadid()]
  @views GradZ = TCacheC1[Threads.threadid()]
  @views FluxZ = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz 
    @. DXpC = 0.0
    @views DerivativeX!(DXpC,pC[:,:,iz],D)
    @. DYpC = 0.0
    @views DerivativeY!(DYpC,pC[:,:,iz],D)
    @views @. FuC[:,:,iz] -= RhoC[:,:,iz] * 
      ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DXpC + 
      (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * DYpC)
    @views @. FvC[:,:,iz] -=  RhoC[:,:,iz] *
      ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DXpC + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DYpC)
    @views @. Fw[:,:,iz] -= RhoC[:,:,iz] * (dXdxI[:,:,1,iz,1,3] * DXpC + dXdxI[:,:,1,iz,2,3] * DYpC)
    @views @. Fw[:,:,iz+1] -= RhoC[:,:,iz] * (dXdxI[:,:,2,iz,1,3] * DXpC + dXdxI[:,:,2,iz,2,3] * DYpC)
  end  
  @inbounds for iz = 2 : Nz 
    @views @. GradZ = 0.5 * (pC[:,:,iz] - pC[:,:,iz-1]) * RhoC[:,:,iz] 
    @views @. FluxZ = GradZ * dXdxI[:,:,1,iz,3,1]
    @views @. FuC[:,:,iz] -= FluxZ 
    @views @. FluxZ = GradZ * dXdxI[:,:,1,iz,3,2]
    @views @. FvC[:,:,iz] -= FluxZ 
    @views @. FluxZ = GradZ * dXdxI[:,:,1,iz,3,3]
    @views @. Fw[:,:,iz] -= FluxZ 
  end    
  @inbounds for iz = 1 : Nz - 1 
    @views @. GradZ = 0.5 * (pC[:,:,iz+1] - pC[:,:,iz]) * RhoC[:,:,iz]  
    @views @. FluxZ = GradZ * dXdxI[:,:,2,iz,3,1]
    @views @. FuC[:,:,iz] -= FluxZ 
    @views @. FluxZ = GradZ * dXdxI[:,:,2,iz,3,2]
    @views @. FvC[:,:,iz] -= FluxZ 
    @views @. FluxZ =  GradZ * dXdxI[:,:,2,iz,3,3]
    @views @. Fw[:,:,iz+1] -= FluxZ 
  end    
  if Nz  == 2
    @views @. GradZ = 0.5 * (pC[:,:,1+1] - pC[:,:,1]) * RhoC[:,:,1]
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,1]
    @views @. FuC[:,:,1] -= FluxZ
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,2]
    @views @. FvC[:,:,1] -= FluxZ
    @views @. FluxZ =  GradZ * dXdxI[:,:,1,1,3,3]
    @views @. Fw[:,:,1] -= FluxZ
  elseif Nz > 2
    for i = 1 : OrdPoly + 1
      for j = 1 : OrdPoly + 1
        @views GradZ[i,j] = (BoundaryDP(pC[i,j,1],pC[i,j,2],pC[i,j,3],
         J[i,j,:,1],J[i,j,:,2],J[i,j,:,3]) - 0.5 * (pC[i,j,1+1] - pC[i,j,1])) * RhoC[i,j,1]
      end
    end
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,1]
    @views @. FuC[:,:,1] -= FluxZ
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,2]
    @views @. FvC[:,:,1] -= FluxZ
    @views @. FluxZ =  GradZ * dXdxI[:,:,1,1,3,3]
    @views @. Fw[:,:,1] -= FluxZ
  end
end 

function GradColumn!(FuC,FvC,Fw,pC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2 = ThreadCache
  Nz = size(FuC,3)
  OrdPoly = Fe.OrdPoly
  D = Fe.DS

  @views DXpC = TCacheC1[Threads.threadid()]
  @views DYpC = TCacheC2[Threads.threadid()]
  @views GradZ = TCacheC1[Threads.threadid()]
  @views FluxZ = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz 
    @. DXpC = 0.0
    @views DerivativeX!(DXpC,pC[:,:,iz],D)
    @. DYpC = 0.0
    @views DerivativeY!(DYpC,pC[:,:,iz],D)
    @views @. FuC[:,:,iz] -=  
      ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DXpC + 
      (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * DYpC)
    @views @. FvC[:,:,iz] -=   
      ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DXpC + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DYpC)
    @views @. Fw[:,:,iz] -= (dXdxI[:,:,1,iz,1,3] * DXpC + dXdxI[:,:,1,iz,2,3] * DYpC)
    @views @. Fw[:,:,iz+1] -= (dXdxI[:,:,2,iz,1,3] * DXpC + dXdxI[:,:,2,iz,2,3] * DYpC)
  end  
  @inbounds for iz = 2 : Nz 
    @views @. GradZ = 0.5 * (pC[:,:,iz] - pC[:,:,iz-1])
    @views @. FluxZ = GradZ * dXdxI[:,:,1,iz,3,1]
    @views @. FuC[:,:,iz] -= FluxZ 
    @views @. FluxZ = GradZ * dXdxI[:,:,1,iz,3,2]
    @views @. FvC[:,:,iz] -= FluxZ 
    @views @. FluxZ = GradZ * dXdxI[:,:,1,iz,3,3]
    @views @. Fw[:,:,iz] -= FluxZ 
  end    
  @inbounds for iz = 1 : Nz - 1 
    @views @. GradZ = 0.5 * (pC[:,:,iz+1] - pC[:,:,iz])
    @views @. FluxZ = GradZ * dXdxI[:,:,2,iz,3,1]
    @views @. FuC[:,:,iz] -= FluxZ 
    @views @. FluxZ = GradZ * dXdxI[:,:,2,iz,3,2]
    @views @. FvC[:,:,iz] -= FluxZ 
    @views @. FluxZ =  GradZ * dXdxI[:,:,2,iz,3,3]
    @views @. Fw[:,:,iz+1] -= FluxZ 
  end    
  if Nz  == 2
    @views @. GradZ = 0.5 * (pC[:,:,1+1] - pC[:,:,1])
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,1]
    @views @. FuC[:,:,1] -= FluxZ
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,2]
    @views @. FvC[:,:,1] -= FluxZ 
    @views @. FluxZ =  GradZ * dXdxI[:,:,1,1,3,3]
    @views @. Fw[:,:,1] -= FluxZ 
  elseif Nz > 2
    for i = 1 : OrdPoly + 1 
      for j = 1 : OrdPoly + 1 
        @views GradZ[i,j] = BoundaryDP(pC[i,j,1],pC[i,j,2],pC[i,j,3],
          J[i,j,:,1],J[i,j,:,2],J[i,j,:,3]) - 0.5 * (pC[i,j,1+1] - pC[i,j,1])  
      end
    end  
#   @views @. GradZ = (pC[:,:,2] - pC[:,:,1]) - 0.5 * (pC[:,:,3] - pC[:,:,2])  
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,1]
    @views @. FuC[:,:,1] -= FluxZ
    @views @. FluxZ = GradZ * dXdxI[:,:,1,1,3,2]
    @views @. FvC[:,:,1] -= FluxZ 
    @views @. FluxZ =  GradZ * dXdxI[:,:,1,1,3,3]
    @views @. Fw[:,:,1] -= FluxZ 
  end 
end 

function DivRhoColumn!(FRhoC,uC,vC,w,RhoC,Fe,dXdxI,ThreadCache)
    @unpack TCacheC1, TCacheC2, TCacheCol1 = ThreadCache
  Nz = size(uC,3)
  OrdPoly = Fe.OrdPoly
  D = Fe.DS

  tempZ = TCacheCol1[Threads.threadid()]
  temp = TCacheC1[Threads.threadid()]
  FluxZ = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,1,iz] = -RhoC[:,:,iz] * (dXdxI[:,:,1,iz,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,1,iz,3,2] * vC[:,:,iz] + dXdxI[:,:,1,iz,3,3] * w[:,:,iz])
    @views @. tempZ[:,:,2,iz] = -RhoC[:,:,iz] * (dXdxI[:,:,2,iz,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,2,iz,3,2] * vC[:,:,iz] + dXdxI[:,:,2,iz,3,3] * w[:,:,iz+1])
  end    
  @inbounds for iz = 1 :Nz  
    @views @. temp = -RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * vC[:,:,iz] + 
      dXdxI[:,:,1,iz,1,3] * w[:,:,iz] + dXdxI[:,:,2,iz,1,3] * w[:,:,iz+1])
    @views DerivativeX!(FRhoC[:,:,iz],temp,D)
    @views @. temp = -RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * uC[:,:,iz] +
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * vC[:,:,iz] + 
      dXdxI[:,:,1,iz,2,3] * w[:,:,iz] + dXdxI[:,:,2,iz,2,3] * w[:,:,iz+1])
    @views DerivativeY!(FRhoC[:,:,iz],temp,D)

    @views @. FRhoC[:,:,iz] += (tempZ[:,:,2,iz] - tempZ[:,:,1,iz])
  end  
  @inbounds for iz = 2 :Nz  
    @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,2,iz-1]) 
    @views @. FRhoC[:,:,iz] += FluxZ
  end  
  @inbounds for iz = 1 :Nz - 1  
    @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] - tempZ[:,:,2,iz]) 
    @views @. FRhoC[:,:,iz] += FluxZ
  end 
end 

function DivRhoTrColumn!(FRhoTrC,uC,vC,w,RhoTrC,Fe,dXdxI,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheCol1 = ThreadCache
  Nz = size(uC,3)
  OrdPoly = Fe.OrdPoly
  D = Fe.DS

  tempZ = TCacheCol1[Threads.threadid()]
  temp = TCacheC1[Threads.threadid()]
  FluxZ = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,1,iz] = -RhoTrC[:,:,iz] * (dXdxI[:,:,1,iz,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,1,iz,3,2] * vC[:,:,iz] + dXdxI[:,:,1,iz,3,3] * w[:,:,iz])
    @views @. tempZ[:,:,2,iz] = -RhoTrC[:,:,iz] * (dXdxI[:,:,2,iz,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,2,iz,3,2] * vC[:,:,iz] + dXdxI[:,:,2,iz,3,3] * w[:,:,iz+1])
  end    
  @inbounds for iz = 1 : Nz  
    @views @. temp = -RhoTrC[:,:,iz] * ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * vC[:,:,iz] + 
      dXdxI[:,:,1,iz,1,3] * w[:,:,iz] + dXdxI[:,:,2,iz,1,3] * w[:,:,iz+1])
    @views DerivativeX!(FRhoTrC[:,:,iz],temp,D)
    @views @. temp = -RhoTrC[:,:,iz] * ((dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * uC[:,:,iz] +
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * vC[:,:,iz] + 
      dXdxI[:,:,1,iz,2,3] * w[:,:,iz] + dXdxI[:,:,2,iz,2,3] * w[:,:,iz+1])
    @views DerivativeY!(FRhoTrC[:,:,iz],temp,D)

    @views @. FRhoTrC[:,:,iz] += (tempZ[:,:,2,iz] - tempZ[:,:,1,iz])
  end    
  @inbounds for iz = 2 : Nz  
    @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,2,iz-1]) 
    @views @. FRhoTrC[:,:,iz] += FluxZ
  end  
  @inbounds for iz = 1 : Nz - 1  
    @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] - tempZ[:,:,2,iz]) 
    @views @. FRhoTrC[:,:,iz] += FluxZ
  end 
end 

function DivRhoGrad!(F,cC,RhoC,Fe,dXdxI,J,ThreadCache,Koeff)
  @unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = ThreadCache
  Nz = size(F,3)
  D = Fe.DS
  DW = Fe.DW

  DxcC = TCacheC1[Threads.threadid()]
  DycC = TCacheC2[Threads.threadid()]
  GradDx = TCacheC3[Threads.threadid()]
  GradDy = TCacheC4[Threads.threadid()]
  temp = TCacheC1[Threads.threadid()]
  Div = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. DxcC = 0.0
    @. DycC = 0.0
    @views DerivativeX!(DxcC,cC[:,:,iz],D)
    @views DerivativeY!(DycC,cC[:,:,iz],D)

    @views @. GradDx = RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DxcC + 
      (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * DycC)
    @views @. GradDy = RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DxcC + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DycC)

    @. Div = 0.0
    @views @. temp = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * GradDx + 
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * GradDy 
    DerivativeX!(Div,temp,DW)
    @views @. temp = (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * GradDx + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * GradDy
    DerivativeY!(Div,temp,DW)

    @views @. F[:,:,iz] -= Koeff * Div / (J[:,:,1,iz] + J[:,:,2,iz])
  end    
end    

function DivRhoGrad!(F,cC,RhoC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = ThreadCache
  Nz = size(F,3)
  D = Fe.DS
  DW = Fe.DW

  DxcC = TCacheC1[Threads.threadid()]
  DycC = TCacheC2[Threads.threadid()]
  GradDx = TCacheC3[Threads.threadid()]
  GradDy = TCacheC4[Threads.threadid()]
  tempC = TCacheC3[Threads.threadid()]
  temp = TCacheC1[Threads.threadid()]
  Div = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. DxcC = 0.0
    @. DycC = 0.0
    @views @. tempC = cC[:,:,iz] / RhoC[:,:,iz]
    DerivativeX!(DxcC,tempC,D)
    DerivativeY!(DycC,tempC,D)

    @views @. GradDx = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DxcC + 
      (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * DycC
    @views @. GradDy = (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DxcC + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DycC

    @. Div = 0.0
    @views @. temp = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * GradDx + 
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * GradDy
    DerivativeX!(Div,temp,DW)
    @views @. temp = (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * GradDx + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * GradDy
    DerivativeY!(Div,temp,DW)

    @views @. F[:,:,iz] = Div / (J[:,:,1,iz] + J[:,:,2,iz])
  end    
end    

function GradDiv!(FuC,FvC,uC,vC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache
  Nz = size(FuC,3)
  D = Fe.DS
  DW = Fe.DW

  Div = TCacheC1[Threads.threadid()]
  temp = TCacheC2[Threads.threadid()]
  DxDiv  = TCacheC2[Threads.threadid()]
  DyDiv  = TCacheC3[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. Div = 0.0  
    @views @. temp = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * vC[:,:,iz]
    DerivativeX!(Div,temp,D)
    @views @. temp = (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * uC[:,:,iz] + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * vC[:,:,iz]
    DerivativeY!(Div,temp,D)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DW)
    DerivativeY!(DyDiv,Div,DW)

    @views @. FuC[:,:,iz] = ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DxDiv +
      (dXdxI[:,:,1,iz,2,1] +dXdxI[:,:,2,iz,2,1]) * DyDiv) / (J[:,:,1,iz] + J[:,:,2,iz])
    @views @. FvC[:,:,iz] = ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DxDiv +
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DyDiv) / (J[:,:,1,iz] + J[:,:,2,iz])
  end  
end  

function GradDiv!(FuC,FvC,uC,vC,RhoC,Fe,dXdxI,J,ThreadCache,Koeff)
  @unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache
  Nz = size(FuC,3)
  D = Fe.DS
  DW = Fe.DW

  Div = TCacheC1[Threads.threadid()]
  temp = TCacheC2[Threads.threadid()]
  DxDiv  = TCacheC2[Threads.threadid()]
  DyDiv  = TCacheC3[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. Div = 0.0  
    @views @. temp = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * vC[:,:,iz]
    DerivativeX!(Div,temp,D)
    @views @. temp = (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * uC[:,:,iz] + 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * vC[:,:,iz]
    DerivativeY!(Div,temp,D)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DW)
    DerivativeY!(DyDiv,Div,DW)

    @views @. FuC[:,:,iz] -= Koeff * RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DxDiv +
      (dXdxI[:,:,1,iz,2,1] +dXdxI[:,:,2,iz,2,1]) * DyDiv) / (J[:,:,1,iz] + J[:,:,2,iz])
    @views @. FvC[:,:,iz] -= Koeff * RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DxDiv +
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DyDiv) / (J[:,:,1,iz] + J[:,:,2,iz])
  end  
end  

function RotCurl!(FuC,FvC,uC,vC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache
  Nz = size(FuC,3)
  D = Fe.DS
  DW = Fe.DW

  W = TCacheC1[Threads.threadid()]
  temp = TCacheC2[Threads.threadid()]
  DxW = TCacheC2[Threads.threadid()]
  DyW = TCacheC3[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. W = 0.0
    @views @. temp = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * uC[:,:,iz]
    DerivativeX!(W,temp,D)
    @views @. temp = (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * uC[:,:,iz]
    DerivativeY!(W,temp,D)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DW)
    DerivativeY!(DyW,W,DW)
    @views @. FvC[:,:,iz] = (-(dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DxW -
      (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * DyW) / (J[:,:,1,iz] + J[:,:,2,iz])
    @views @. FuC[:,:,iz] = ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DxW +
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DyW) / (J[:,:,1,iz] + J[:,:,2,iz])
  end   
end   

function RotCurl!(FuC,FvC,uC,vC,RhoC,Fe,dXdxI,J,ThreadCache,Koeff)
  @unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache
  Nz = size(FuC,3)
  D = Fe.DS
  DW = Fe.DW

  W = TCacheC1[Threads.threadid()]
  temp = TCacheC2[Threads.threadid()]
  DxW = TCacheC2[Threads.threadid()]
  DyW = TCacheC3[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. W = 0.0
    @views @. temp = (dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * uC[:,:,iz]
    DerivativeX!(W,temp,D)
    @views @. temp = (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * uC[:,:,iz]
    DerivativeY!(W,temp,D)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DW)
    DerivativeY!(DyW,W,DW)
    @views @. FvC[:,:,iz] -= Koeff * RhoC[:,:,iz] * (-(dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * DxW -
      (dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * DyW) / (J[:,:,1,iz] + J[:,:,2,iz])
    @views @. FuC[:,:,iz] -= Koeff * RhoC[:,:,iz] * ((dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * DxW +
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * DyW) / (J[:,:,1,iz] + J[:,:,2,iz])
  end   
end   

function Buoyancy!(FwF,RhoC,J,Phys)

  @views @. FwF[:,:,2:end-1] -= Phys.Grav *
    (RhoC[:,:,1:end-1] * J[:,:,2,1:end-1] +
    RhoC[:,:,2:end] * J[:,:,1,2:end])

end

function Rot!(Rot,uC,vC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheC3 = ThreadCache

  Nz = size(Rot,3)
  D = Fe.DS

  W = TCacheC1[Threads.threadid()]
  temp = TCacheC2[Threads.threadid()]

  @inbounds for iz = 1 : Nz
    @. W = 0.0
    @views @. temp = 0.5 * ((dXdxI[:,:,1,iz,1,1] + dXdxI[:,:,2,iz,1,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,1,iz,1,2] + dXdxI[:,:,2,iz,1,2]) * uC[:,:,iz])
    DerivativeX!(W,temp,D)
    @views @. temp = 0.5 * ((dXdxI[:,:,1,iz,2,1] + dXdxI[:,:,2,iz,2,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,1,iz,2,2] + dXdxI[:,:,2,iz,2,2]) * uC[:,:,iz])
    DerivativeY!(W,temp,D)
    @views @. Rot[:,:,iz] = 2.0 * W / (J[:,:,1,iz] + J[:,:,2,iz])
  end  
end

function BoundaryW!(wCG,v1CG,v2CG,CG,dXdxI)
  @views @. wCG = -(dXdxI[:,:,3,1] * v1CG +
    dXdxI[:,:,3,2] * v2CG)/
    dXdxI[:,:,3,3]
end

function BoundaryDP(p1,p2,p3,J1,J2,J3)

  h1 = 0.5 * (J1[1] + J1[2]);
  h2 = 0.5 * (J2[1] + J2[2]);
  h3 = 0.5 * (J3[1] + J3[2]);
  f1 = -h1 * (3*h1 + 4*h2 + 2*h3)/((h1 + h2)*(h1 + h2 + h3)) ;
  f2 = h1 * (h1^2 + 6*h1*h2 + 3*h1*h3 + 6*h2^2 + 6*h2*h3 + 2*h3^2)/ 
   ((h1 + h2)*(h2 + h3)*(h1 + h2 + h3)) 
  f3 = -h1 * (h1 + 2*h2)/((h2 + h3)*(h1 + h2 + h3)) 
  DpB = (f1 * p1 + f2 * p2 + f3 * p3);

end  


