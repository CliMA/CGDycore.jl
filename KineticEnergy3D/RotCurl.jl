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
    @views @. temp = cC[iz,:,:] / RhoC[iz,:,:]
    DerivativeX!(DxcC,temp,DX)
    DerivativeY!(DycC,temp,DY)

    @views @. GradDx = 0.5 * RhoC[iz,:,:,] * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DxcC + 
      (dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DycC)
    @views @. GradDy = 0.5 * RhoC[iz,:,:,] * ((dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * DxcC + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DycC)

    @. Div = 0.0
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * GradDx + 
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * GradDy)
    DerivativeX!(Div,temp,DXW)
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * GradDx + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * GradDy)
    DerivativeY!(Div,temp,DYW)

    @views @. F[iz,:,:] -= 2.0 * Koeff * Div / (J[iz,:,:,1] + J[iz,:,:,2])
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
    @views @. temp = cC[iz,:,:] / RhoC[iz,:,:]
    DerivativeX!(DxcC,temp,DX)
    DerivativeY!(DycC,temp,DY)

    @views @. GradDx = 0.5 *((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DxcC + 
      (dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DycC)
    @views @. GradDy = 0.5 *((dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * DxcC + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DycC)

    @. Div = 0.0
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * GradDx + 
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * GradDy)
    DerivativeX!(Div,temp,DXW)
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * GradDx + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * GradDy)
    DerivativeY!(Div,temp,DYW)

    @views @. F[iz,:,:] = 2.0 * Div / (J[iz,:,:,1] + J[iz,:,:,2])
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
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,2,1,1,1]) * uC[iz,:,:] +
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * vC[iz,:,:])
    DerivativeX!(Div,temp,DX)
    @views @. temp = 0.5 *((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * uC[iz,:,:] + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * vC[iz,:,:])
    DerivativeY!(Div,temp,DY)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DXW)
    DerivativeY!(DyDiv,Div,DYW)

    @views @. FuC[iz,:,:] += ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DxDiv +
      (dXdxI[iz,:,:,1,2,1] +dXdxI[iz,:,:,2,2,1]) * DyDiv) / (J[iz,:,:,1] + J[iz,:,:,2])
    @views @. FvC[iz,:,:] += ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DxDiv +
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DyDiv) / (J[iz,:,:,1] + J[iz,:,:,2])
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
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,2,1,1,1]) * uC[iz,:,:] +
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * vC[iz,:,:])
    DerivativeX!(Div,temp,DX)
    @views @. temp = 0.5 *((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * uC[iz,:,:] + 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * vC[iz,:,:])
    DerivativeY!(Div,temp,DY)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DXW)
    DerivativeY!(DyDiv,Div,DYW)

    @views @. FuC[iz,:,:] -= Koeff * RhoC[iz,:,:] * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DxDiv +
      (dXdxI[iz,:,:,1,2,1] +dXdxI[iz,:,:,2,2,1]) * DyDiv) / (J[iz,:,:,1] + J[iz,:,:,2])
    @views @. FvC[iz,:,:] -= Koeff * RhoC[iz,:,:] * ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DxDiv +
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DyDiv) / (J[iz,:,:,1] + J[iz,:,:,2])
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
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * vC[iz,:,:] - 
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * uC[iz,:,:])
    DerivativeX!(W,temp,DX)
    @views @. temp = 0.5 * (dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * vC[iz,:,:] - 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * uC[iz,:,:]
    DerivativeY!(W,temp,DY)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DXW)
    DerivativeY!(DyW,W,DYW)
    @views @. FvC[iz,:,:] += (-(dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DxW -
      (dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DyW) / (J[iz,:,:,1] + J[iz,:,:,2])
    @views @. FuC[iz,:,:] += ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DxW +
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DyW) / (J[iz,:,:,1] + J[iz,:,:,2])
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
    @views @. temp = 0.5 * ((dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * vC[iz,:,:] - 
      (dXdxI[iz,:,:,1,1,2] + dXdxI[iz,:,:,2,1,2]) * uC[iz,:,:])
    DerivativeX!(W,temp,DX)
    @views @. temp = 0.5 * (dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * vC[iz,:,:] - 
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * uC[iz,:,:]
    DerivativeY!(W,temp,DY)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DXW)
    DerivativeY!(DyW,W,DYW)
    @views @. FvC[iz,:,:] -= Koeff * RhoC[iz,:,:] * (-(dXdxI[iz,:,:,1,1,1] + dXdxI[iz,:,:,2,1,1]) * DxW -
      (dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DyW) / (J[iz,:,:,1] + J[iz,:,:,2])
    @views @. FuC[iz,:,:] -= Koeff * RhoC[iz,:,:] * ((dXdxI[iz,:,:,1,2,1] + dXdxI[iz,:,:,2,2,1]) * DxW +
      (dXdxI[iz,:,:,1,2,2] + dXdxI[iz,:,:,2,2,2]) * DyW) / (J[iz,:,:,1] + J[iz,:,:,2])
  end   
end   
