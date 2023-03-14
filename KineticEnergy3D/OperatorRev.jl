function CurlColumn!(FuC,FvC,Fw,uC,vC,w,RhoC,Fe,dXdxI,Cache)
  Nz = size(FuC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempuZ = Cache.Column[:,:,1:2,:,1]
  @views tempvZ = Cache.Column[:,:,1:2,:,2]
  @views tempwZ = Cache.Column[:,:,1:2,:,3]
  @views U = Cache.Block[:,:,:,1]
  @views V = Cache.Block[:,:,:,2]
  @views W = Cache.Block[:,:,:,3]
  @views temp = Cache.Block[:,:,1:2,4]

  @views FluxUZ = Cache.BlockXY[:,:,1]
  @views FluxVZ = Cache.BlockXY[:,:,2]
  @views FluxWZ = Cache.BlockXY[:,:,3]
  OPz = OrdPolyZ + 1
  @inbounds for iz = 1 : Nz 
    @views @. tempuZ[:,:,1,iz] = -dXdxI[:,:,iz,1,3,3] * vC[:,:,iz] + dXdxI[:,:,iz,1,3,2] * w[:,:,iz]
    @views @. tempuZ[:,:,2,iz] = -dXdxI[:,:,iz,2,3,3] * vC[:,:,iz] + dXdxI[:,:,iz,2,3,2] * w[iz+1,:,:]
    @views @. tempvZ[:,:,1,iz] = dXdxI[:,:,iz,1,3,3] * uC[:,:,iz] - dXdxI[:,:,iz,1,3,1] * w[:,:,iz]
    @views @. tempvZ[:,:,2,iz] = dXdxI[:,:,iz,2,3,3] * uC[:,:,iz] - dXdxI[:,:,iz,2,3,1] * w[iz+1,:,:]
    @views @. tempwZ[:,:,1,iz] = dXdxI[:,:,iz,1,3,1] * vC[:,:,iz] - dXdxI[:,:,iz,1,3,2] * uC[:,:,iz]
    @views @. tempwZ[:,:,2,iz] = dXdxI[:,:,iz,2,3,1] * vC[:,:,iz] - dXdxI[:,:,iz,2,3,2] * uC[:,:,iz]
  end  
  @inbounds for iz = 1 : Nz 
    if iz > 1
      @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz] - tempuZ[:,:,2,iz-1])
      @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz] - tempvZ[:,:,2,iz-1])
      @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz] - tempwZ[:,:,2,iz-1])
      @views @. FuC[:,:,iz] += RhoC[:,:,iz] * (vC[:,:,iz,1] * FluxWZ - w[:,:,iz,1] * FluxVZ)
      @views @. FvC[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz,1] * FluxWZ - w[:,:,iz,1] * FluxUZ)
      @views @. Fw[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz,1] * FluxVZ - vC[:,:,iz,1] * FluxUZ)
    end
    if iz < Nz
      @views @. FluxUZ = 0.5 * (tempuZ[:,:,1,iz+1] - tempuZ[:,:,2,iz])
      @views @. FluxVZ = 0.5 * (tempvZ[:,:,1,iz+1] - tempvZ[:,:,2,iz])
      @views @. FluxWZ = 0.5 * (tempwZ[:,:,1,iz+1] - tempwZ[:,:,2,iz])
      @views @. FuC[:,:,iz] += RhoC[:,:,iz] * (vC[:,:,iz,] * FluxWZ - w[iz+1,:,:] * FluxVZ)
      @views @. FvC[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz,] * FluxWZ - w[iz+1,:,:] * FluxUZ)
      @views @. Fw[iz+1,:,:] += RhoC[:,:,iz,] * (uC[:,:,iz] * FluxVZ - vC[:,:,iz] * FluxUZ)
    end    
  
    @. U = 0.0
    @. V = 0.0
    @. W = 0.0
    @views @. temp[:,:,1] = -dXdxI[:,:,iz,1,1,3] * vC[:,:,iz] + dXdxI[:,:,iz,1,1,2] * w[:,:,iz]
    @views @. temp[:,:,2] = -dXdxI[:,:,iz,2,1,3] * vC[:,:,iz] + dXdxI[:,:,iz,2,1,2] * w[iz+1,:,:]
    DerivativeX!(U,temp,DX)
    @views @. temp[:,:,1] = dXdxI[:,:,iz,1,1,3] * uC[:,:,iz] - dXdxI[:,:,iz,1,1,1] * w[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,iz,2,1,3] * uC[:,:,iz] - dXdxI[:,:,iz,2,1,1] * w[iz+1,:,:]
    DerivativeX!(V,temp,DX)
    @views @. temp[:,:,1] = dXdxI[:,:,iz,1,1,1] * vC[:,:,iz] - dXdxI[:,:,iz,1,1,2] * uC[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,iz,2,1,1] * vC[:,:,iz] - dXdxI[:,:,iz,2,1,2] * uC[:,:,iz]
    DerivativeX!(W,temp,DX)
    @views @. temp[:,:,1] = -dXdxI[:,:,iz,1,2,3] * vC[:,:,iz] + dXdxI[:,:,iz,1,2,2] * w[:,:,iz]
    @views @. temp[:,:,2] = -dXdxI[:,:,iz,2,2,3] * vC[:,:,iz] + dXdxI[:,:,iz,2,2,2] * w[iz+1,:,:]
    DerivativeY!(U,temp,DY)
    @views @. temp[:,:,1] = dXdxI[:,:,iz,1,2,3] * uC[:,:,iz] - dXdxI[:,:,iz,1,2,1] * w[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,iz,2,2,3] * uC[:,:,iz] - dXdxI[:,:,iz,2,2,1] * w[iz+1,:,:]
    DerivativeY!(V,temp,DY)
    @views @. temp[:,:,1] = dXdxI[:,:,iz,1,2,1] * vC[:,:,iz] - dXdxI[:,:,iz,1,2,2] * uC[:,:,iz]
    @views @. temp[:,:,2] = dXdxI[:,:,iz,2,2,1] * vC[:,:,iz] - dXdxI[:,:,iz,2,2,2] * uC[:,:,iz]
    DerivativeY!(W,temp,DY)
    @views DerivativeZ!(U,tempuZ[:,:,:,iz],DZ)
    @views DerivativeZ!(V,tempvZ[:,:,:,iz],DZ)
    @views DerivativeZ!(W,tempwZ[:,:,:,iz],DZ)
    @views @. FuC[:,:,iz] += RhoC[:,:,iz] * (vC[:,:,iz] * W[:,:,1] - w[:,:,iz] * V[:,:,1] +
      vC[:,:,iz] * W[:,:,2] - w[iz+1,:,:] * V[:,:,2]) 
    @views @. FvC[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz] * W[:,:,1] - w[:,:,iz,:] * U[:,:,1] +
      uC[:,:,iz] * W[:,:,2] - w[iz+1,:,:,:] * U[:,:,2]) 
    @views @. Fw[:,:,iz] += RhoC[:,:,iz] * (uC[:,:,iz] * V[:,:,1] - vC[:,:,iz] * U[:,:,1])
    @views @. Fw[iz+1,:,:] += RhoC[:,:,iz] * (uC[:,:,iz] * V[:,:,2] - vC[:,:,iz] * U[:,:,2])
      
  end 
end

function RhoGradColumn!(FuC,FvC,Fw,pC,RhoC,Fe,dXdxI,Cache)
  Nz = size(pC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views DXpC = Cache.Block[:,:,1,1]
  @views DYpC = Cache.Block[:,:,1,2]
  @views DZpC = Cache.Block[:,:,1,3]
  @views GradZ = Cache.BlockXY[:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,2]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz 
    @. DXpC = 0.0
    @views DerivativeX!(DXpC,pC[:,:,iz],DX)
    @. DYpC = 0.0
    @views DerivativeY!(DYpC,pC[:,:,iz],DY)
#   @. DZpC = 0.0
#   @views DerivativeZ!(DZpC,pC[:,:,iz],DZ)
    @views @. FuC[:,:,iz] -= RhoC[:,:,iz] * 
      ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DXpC + 
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,2,1]) * DYpC)
    @views @. FvC[:,:,iz] -=  RhoC[:,:,iz] *
      ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DXpC + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DYpC)
    @views @. Fw[:,:,iz] -= RhoC[:,:,iz] * (dXdxI[:,:,iz,1,1,3] * DXpC + dXdxI[:,:,iz,1,2,3] * DYpC)
    @views @. Fw[iz+1,:,:] -= RhoC[:,:,iz] * (dXdxI[:,:,iz,2,1,3] * DXpC + dXdxI[:,:,iz,2,2,3] * DYpC)
    if iz > 1
      @views @. GradZ = 0.5 * (pC[:,:,iz] - pC[iz-1,:,:]) * RhoC[:,:,iz] 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,1,3,1]
      @views @. FuC[:,:,iz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,1,3,2]
      @views @. FvC[:,:,iz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,1,3,3]
      @views @. Fw[:,:,iz] -= FluxZ 
    end    
    if iz < Nz 
      @views @. GradZ = 0.5 * (pC[iz+1,:,:] - pC[:,:,iz]) * RhoC[:,:,iz]  
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,2,3,1]
      @views @. FuC[:,:,iz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,2,3,2]
      @views @. FvC[:,:,iz] -= FluxZ 
      @views @. FluxZ =  GradZ * dXdxI[:,:,iz,2,3,3]
      @views @. Fw[iz+1,:,:] -= FluxZ 
    end    
  end 
  @views @. GradZ = 0.5 * (pC[1+1,:,:] - pC[1,:,:]) * RhoC[1,:,:]  
  @views @. FluxZ = GradZ * dXdxI[1,:,:,1,3,1]
  @views @. FuC[1,:,:] -= FluxZ 
  @views @. FluxZ = GradZ * dXdxI[1,:,:,1,3,2]
  @views @. FvC[1,:,:] -= FluxZ 
  @views @. FluxZ =  GradZ * dXdxI[1,:,:,1,3,3]
  @views @. Fw[1,:,:] -= FluxZ 
end 

function GradColumn!(FuC,FvC,Fw,pC,Fe,dXdxI,Cache)
  Nz = size(pC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views DXpC = Cache.Block[:,:,1,1]
  @views DYpC = Cache.Block[:,:,1,2]
  @views DZpC = Cache.Block[:,:,1,3]
  @views GradZ = Cache.BlockXY[:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,2]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz 
    @. DXpC = 0.0
    @views DerivativeX!(DXpC,pC[:,:,iz],DX)
    @. DYpC = 0.0
    @views DerivativeY!(DYpC,pC[:,:,iz],DY)
    @views @. FuC[:,:,iz] -=  
      ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DXpC + 
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,2,1]) * DYpC)
    @views @. FvC[:,:,iz] -=   
      ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DXpC + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DYpC)
    @views @. Fw[:,:,iz] -= (dXdxI[:,:,iz,1,1,3] * DXpC + dXdxI[:,:,iz,1,2,3] * DYpC)
    @views @. Fw[iz+1,:,:] -= (dXdxI[:,:,iz,2,1,3] * DXpC + dXdxI[:,:,iz,2,2,3] * DYpC)
    if iz > 1
      @views @. GradZ = 0.5 * (pC[:,:,iz] - pC[iz-1,:,:])
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,1,3,1]
      @views @. FuC[:,:,iz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,1,3,2]
      @views @. FvC[:,:,iz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,1,3,3]
      @views @. Fw[:,:,iz] -= FluxZ 
    end    
    if iz < Nz 
      @views @. GradZ = 0.5 * (pC[iz+1,:,:] - pC[:,:,iz])
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,2,3,1]
      @views @. FuC[:,:,iz] -= FluxZ 
      @views @. FluxZ = GradZ * dXdxI[:,:,iz,2,3,2]
      @views @. FvC[:,:,iz] -= FluxZ 
      @views @. FluxZ =  GradZ * dXdxI[:,:,iz,2,3,3]
      @views @. Fw[iz+1,:,:] -= FluxZ 
    end    
  end 
  @views @. GradZ = 0.5 * (pC[1+1,:,:] - pC[1,:,:])
  @views @. FluxZ = GradZ * dXdxI[1,:,:,1,3,1]
  @views @. FuC[1,:,:] -= FluxZ
  @views @. FluxZ = GradZ * dXdxI[1,:,:,1,3,2]
  @views @. FvC[1,:,:] -= FluxZ 
  @views @. FluxZ =  GradZ * dXdxI[1,:,:,1,3,3]
  @views @. Fw[1,:,:] -= FluxZ 
end 

function DivRhoColumn!(FRhoC,uC,vC,w,RhoC,Fe,dXdxI,Cache)
  Nz = size(uC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempZ = Cache.Column[:,:,:,:,1]
  @views temp = Cache.Block[:,:,1,1]
  @views temp1 = Cache.Block[:,:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,1]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,1,iz] = -RhoC[:,:,iz] * (dXdxI[:,:,iz,1,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,iz,1,3,2] * vC[:,:,iz] + dXdxI[:,:,iz,1,3,3] * w[:,:,iz])
    @views @. tempZ[:,:,2,iz] = -RhoC[:,:,iz] * (dXdxI[:,:,iz,2,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,iz,2,3,2] * vC[:,:,iz] + dXdxI[:,:,iz,2,3,3] * w[iz+1,:,:])
  end    
  @inbounds for iz = 1 :Nz  
    @views @. temp = -RhoC[:,:,iz] * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * vC[:,:,iz] + 
      dXdxI[:,:,iz,1,1,3] * w[:,:,iz] + dXdxI[:,:,iz,2,1,3] * w[iz+1,:,:])
    @views DerivativeX!(FRhoC[:,:,iz],temp,DX)
    @views @. temp = -RhoC[:,:,iz] * ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * uC[:,:,iz] +
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * vC[:,:,iz] + 
      dXdxI[:,:,iz,1,2,3] * w[:,:,iz] + dXdxI[:,:,iz,2,2,3] * w[iz+1,:,:])
    @views DerivativeY!(FRhoC[:,:,iz],temp,DY)
    @views @. temp1 = 0.0 
    @views DerivativeZ!(temp1,tempZ[:,:,:,iz],DZ)
    @views @. FRhoC[:,:,iz] += (temp1[:,:,1] + temp1[:,:,2]) 
    # Density
    if iz > 1
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,2,iz-1]) 
      @views @. FRhoC[:,:,iz] += FluxZ
    end  
    if iz < Nz
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] - tempZ[:,:,2,iz]) 
      @views @. FRhoC[:,:,iz] += FluxZ
    end  
  end 
end 

function DivRhoThetaColumn!(FRhoThetaC,uC,vC,w,RhoThetaC,Fe,dXdxI,Cache)
  Nz = size(uC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  DX = Fe.DX
  DY = Fe.DY
  DZ = Fe.DZ

  @views tempZ = Cache.Column[:,:,:,:,1]
  @views temp = Cache.Block[:,:,1,1]
  @views temp1 = Cache.Block[:,:,:,1]
  @views FluxZ = Cache.BlockXY[:,:,1]

  OPz = OrdPolyZ + 1

  @inbounds for iz = 1 : Nz  
    @views @. tempZ[:,:,1,iz] = -RhoThetaC[:,:,iz] * (dXdxI[:,:,iz,1,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,iz,1,3,2] * vC[:,:,iz] + dXdxI[:,:,iz,1,3,3] * w[:,:,iz])
    @views @. tempZ[:,:,2,iz] = -RhoThetaC[:,:,iz] * (dXdxI[:,:,iz,2,3,1] * uC[:,:,iz] +
      +dXdxI[:,:,iz,2,3,2] * vC[:,:,iz] + dXdxI[:,:,iz,2,3,3] * w[iz+1,:,:])
  end    
  @inbounds for iz = 1 :Nz  
    @views @. temp = -RhoThetaC[:,:,iz] * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * vC[:,:,iz] + 
      dXdxI[:,:,iz,1,1,3] * w[:,:,iz] + dXdxI[:,:,iz,2,1,3] * w[iz+1,:,:])
    @views DerivativeX!(FRhoThetaC[:,:,iz],temp,DX)
    @views @. temp = -RhoThetaC[:,:,iz] * ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * uC[:,:,iz] +
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * vC[:,:,iz] + 
      dXdxI[:,:,iz,1,2,3] * w[:,:,iz] + dXdxI[:,:,iz,2,2,3] * w[iz+1,:,:])
    @views DerivativeY!(FRhoThetaC[:,:,iz],temp,DY)
    @views @. temp1 = 0.0 
    @views DerivativeZ!(temp1,tempZ[:,:,:,iz],DZ)
    @views @. FRhoThetaC[:,:,iz] += (temp1[:,:,1] + temp1[:,:,2]) 
    # Density
    if iz > 1
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz] - tempZ[:,:,2,iz-1]) 
      @views @. FRhoThetaC[:,:,iz] += FluxZ
    end  
    if iz < Nz
      @views @. FluxZ = 0.5 * (tempZ[:,:,1,iz+1] - tempZ[:,:,2,iz]) 
      @views @. FRhoThetaC[:,:,iz] += FluxZ
    end  
  end 
end 

function FDivRhoGrad!(F,cC,RhoC,Fe,dXdxI,J,Cache,Koeff)

  Nz = size(F,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views GradDx = Cache.Block[:,:,1,1]
  @views GradDy = Cache.Block[:,:,1,2]
  @views DxcC = Cache.Block[:,:,1,3]
  @views DycC = Cache.Block[:,:,1,4]
  @views temp = Cache.Block[:,:,1,1]
  @views Div = Cache.Block[:,:,1,4]

  @inbounds for iz = 1 : Nz
    @. DxcC = 0.0
    @. DycC = 0.0
    @views @. temp = cC[:,:,iz] / RhoC[:,:,iz]
    DerivativeX!(DxcC,temp,DX)
    DerivativeY!(DycC,temp,DY)

    @views @. GradDx = 0.5 * RhoC[:,:,iz,] * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DxcC + 
      (dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DycC)
    @views @. GradDy = 0.5 * RhoC[:,:,iz,] * ((dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * DxcC + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DycC)

    @. Div = 0.0
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * GradDx + 
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * GradDy)
    DerivativeX!(Div,temp,DXW)
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * GradDx + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * GradDy)
    DerivativeY!(Div,temp,DYW)

    @views @. F[:,:,iz] -= 2.0 * Koeff * Div / (J[:,:,iz,1] + J[:,:,iz,2])
  end    
end    

function FDivRhoGrad!(F,cC,RhoC,Fe,dXdxI,J,Cache)

  Nz = size(F,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views GradDx = Cache.Block[:,:,1,1]
  @views GradDy = Cache.Block[:,:,1,2]
  @views DxcC = Cache.Block[:,:,1,3]
  @views DycC = Cache.Block[:,:,1,4]
  @views temp = Cache.Block[:,:,1,1]
  @views Div = Cache.Block[:,:,1,4]

  @inbounds for iz = 1 : Nz
    @. DxcC = 0.0
    @. DycC = 0.0
    @views @. temp = cC[:,:,iz] / RhoC[:,:,iz]
    DerivativeX!(DxcC,temp,DX)
    DerivativeY!(DycC,temp,DY)

    @views @. GradDx = 0.5 *((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DxcC + 
      (dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DycC)
    @views @. GradDy = 0.5 *((dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * DxcC + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DycC)

    @. Div = 0.0
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * GradDx + 
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * GradDy)
    DerivativeX!(Div,temp,DXW)
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * GradDx + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * GradDy)
    DerivativeY!(Div,temp,DYW)

    @views @. F[:,:,iz] = 2.0 * Div / (J[:,:,iz,1] + J[:,:,iz,2])
  end    
end    

function GradDiv!(FuC,FvC,uC,vC,Fe,dXdxI,J,Cache)
  Nz = size(FuC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views Div = Cache.Block[:,:,1,1]
  @views temp = Cache.Block[:,:,1,2]
  @views DxDiv = Cache.Block[:,:,1,2]
  @views DyDiv = Cache.Block[:,:,1,3]

  @inbounds for iz = 1 : Nz
    @. Div = 0.0  
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,1,1] + dXdxI[iz,:,2,1,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * vC[:,:,iz])
    DerivativeX!(Div,temp,DX)
    @views @. temp = 0.5 *((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * uC[:,:,iz] + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * vC[:,:,iz])
    DerivativeY!(Div,temp,DY)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DXW)
    DerivativeY!(DyDiv,Div,DYW)

    @views @. FuC[:,:,iz] += ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DxDiv +
      (dXdxI[:,:,iz,1,2,1] +dXdxI[:,:,iz,2,2,1]) * DyDiv) / (J[:,:,iz,1] + J[:,:,iz,2])
    @views @. FvC[:,:,iz] += ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DxDiv +
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DyDiv) / (J[:,:,iz,1] + J[:,:,iz,2])
  end  
end  
function GradDiv!(FuC,FvC,uC,vC,RhoC,Fe,dXdxI,J,Cache,Koeff)
  Nz = size(FuC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views Div = Cache.Block[:,:,1,1]
  @views temp = Cache.Block[:,:,1,2]
  @views DxDiv = Cache.Block[:,:,1,2]
  @views DyDiv = Cache.Block[:,:,1,3]

  @inbounds for iz = 1 : Nz
    @. Div = 0.0  
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,1,1] + dXdxI[iz,:,2,1,1,1]) * uC[:,:,iz] +
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * vC[:,:,iz])
    DerivativeX!(Div,temp,DX)
    @views @. temp = 0.5 *((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * uC[:,:,iz] + 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * vC[:,:,iz])
    DerivativeY!(Div,temp,DY)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DXW)
    DerivativeY!(DyDiv,Div,DYW)

    @views @. FuC[:,:,iz] -= Koeff * RhoC[:,:,iz] * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DxDiv +
      (dXdxI[:,:,iz,1,2,1] +dXdxI[:,:,iz,2,2,1]) * DyDiv) / (J[:,:,iz,1] + J[:,:,iz,2])
    @views @. FvC[:,:,iz] -= Koeff * RhoC[:,:,iz] * ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DxDiv +
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DyDiv) / (J[:,:,iz,1] + J[:,:,iz,2])
  end  
end  

function RotCurl!(FuC,FvC,uC,vC,Fe,dXdxI,J,Cache)
  Nz = size(FuC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views W = Cache.Block[:,:,1,1]
  @views temp = Cache.Block[:,:,1,2]
  @views DxW = Cache.Block[:,:,1,2]
  @views DyW = Cache.Block[:,:,1,3]

  @inbounds for iz = 1 : Nz
    @. W = 0.0
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * uC[:,:,iz])
    DerivativeX!(W,temp,DX)
    @views @. temp = 0.5 * (dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * uC[:,:,iz]
    DerivativeY!(W,temp,DY)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DXW)
    DerivativeY!(DyW,W,DYW)
    @views @. FvC[:,:,iz] += (-(dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DxW -
      (dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DyW) / (J[:,:,iz,1] + J[:,:,iz,2])
    @views @. FuC[:,:,iz] += ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DxW +
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DyW) / (J[:,:,iz,1] + J[:,:,iz,2])
  end   
end   

function RotCurl!(FuC,FvC,uC,vC,RhoC,Fe,dXdxI,J,Cache,Koeff)
  Nz = size(FuC,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views W = Cache.Block[:,:,1,1]
  @views temp = Cache.Block[:,:,1,2]
  @views DxW = Cache.Block[:,:,1,2]
  @views DyW = Cache.Block[:,:,1,3]

  @inbounds for iz = 1 : Nz
    @. W = 0.0
    @views @. temp = 0.5 * ((dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,iz,1,1,2] + dXdxI[:,:,iz,2,1,2]) * uC[:,:,iz])
    DerivativeX!(W,temp,DX)
    @views @. temp = 0.5 * (dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * vC[:,:,iz] - 
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * uC[:,:,iz]
    DerivativeY!(W,temp,DY)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DXW)
    DerivativeY!(DyW,W,DYW)
    @views @. FvC[:,:,iz] -= Koeff * RhoC[:,:,iz] * (-(dXdxI[:,:,iz,1,1,1] + dXdxI[:,:,iz,2,1,1]) * DxW -
      (dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DyW) / (J[:,:,iz,1] + J[:,:,iz,2])
    @views @. FuC[:,:,iz] -= Koeff * RhoC[:,:,iz] * ((dXdxI[:,:,iz,1,2,1] + dXdxI[:,:,iz,2,2,1]) * DxW +
      (dXdxI[:,:,iz,1,2,2] + dXdxI[:,:,iz,2,2,2]) * DyW) / (J[:,:,iz,1] + J[:,:,iz,2])
  end   
end   
