function FDivRhoGrad!(F,cF,RhoF,Fe,dXdxI,J,Cache,Koeff)

  Nz = size(F,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views GradDx = Cache.Block[:,:,:,1]
  @views GradDy = Cache.Block[:,:,:,2]
  @views DxcF = Cache.Block[:,:,:,3]
  @views DycF = Cache.Block[:,:,:,4]
  @views temp = Cache.Block[:,:,:,1]
  @views Div = Cache.Block[:,:,:,4]

  @inbounds for iz = 1 : Nz
    @. DxcF = 0.0
    @. DycF = 0.0
    @views DerivativeX!(DxcF,cF[iz,:,:,:],DX)
    @views DerivativeY!(DycF,cF[iz,:,:,:],DY)

    @views @. GradDx = RhoF[iz,:,:,:] * (dXdxI[iz,:,:,:,1,1] * DxcF + dXdxI[iz,:,:,:,2,1] * DycF)
    @views @. GradDy = RhoF[iz,:,:,:] * (dXdxI[iz,:,:,:,1,2] * DxcF + dXdxI[iz,:,:,:,2,2] * DycF)

    @. Div = 0.0
    @views @. temp = dXdxI[iz,:,:,:,1,1] * GradDx + dXdxI[iz,:,:,:,1,2] * GradDy
    DerivativeX!(Div,temp,DXW)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * GradDx + dXdxI[iz,:,:,:,2,2] * GradDy
    DerivativeY!(Div,temp,DYW)

    @views @. F[iz,:,:,:] -= Koeff * Div / J[iz,:,:,:]
  end    
end    
function FDivRhoGrad!(F,cF,RhoF,Fe,dXdxI,J,Cache)

  Nz = size(F,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views GradDx = Cache.Block[:,:,:,1]
  @views GradDy = Cache.Block[:,:,:,2]
  @views DxcF = Cache.Block[:,:,:,3]
  @views DycF = Cache.Block[:,:,:,4]
  @views temp = Cache.Block[:,:,:,1]
  @views Div = Cache.Block[:,:,:,4]

  @inbounds for iz = 1 : Nz
    @. DxcF = 0.0
    @. DycF = 0.0
    @views @. temp = cF[iz,:,:,:] / RhoF[iz,:,:,:]
    DerivativeX!(DxcF,temp,DX)
    DerivativeY!(DycF,temp,DY)

    @views @. GradDx = dXdxI[iz,:,:,:,1,1] * DxcF + dXdxI[iz,:,:,:,2,1] * DycF
    @views @. GradDy = dXdxI[iz,:,:,:,1,2] * DxcF + dXdxI[iz,:,:,:,2,2] * DycF

    @. Div = 0.0
    @views @. temp = dXdxI[iz,:,:,:,1,1] * GradDx + dXdxI[iz,:,:,:,1,2] * GradDy
    DerivativeX!(Div,temp,DXW)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * GradDx + dXdxI[iz,:,:,:,2,2] * GradDy
    DerivativeY!(Div,temp,DYW)

    @views @. F[iz,:,:,:] = Div / J[iz,:,:,:]
  end    
end    

function GradDiv!(FuF,FvF,uF,vF,Fe,dXdxI,J,Cache)
  Nz = size(FuF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views Div = Cache.Block[:,:,:,1]
  @views temp = Cache.Block[:,:,:,2]
  @views DxDiv = Cache.Block[:,:,:,2]
  @views DyDiv = Cache.Block[:,:,:,3]

  @inbounds for iz = 1 : Nz
    @. Div = 0.0  
    @views @. temp = dXdxI[iz,:,:,:,1,1] * uF[iz,:,:,:] + dXdxI[iz,:,:,:,1,2] * vF[iz,:,:,:]
    DerivativeX!(Div,temp,DX)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * uF[iz,:,:,:] + dXdxI[iz,:,:,:,2,2] * vF[iz,:,:,:]
    DerivativeY!(Div,temp,DY)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DXW)
    DerivativeY!(DyDiv,Div,DYW)

    @views @. FuF[iz,:,:,:] += (dXdxI[iz,:,:,:,1,1] * DxDiv +
      dXdxI[iz,:,:,:,2,1] * DyDiv) / J[iz,:,:,:]
    @views @. FvF[iz,:,:,:] += (dXdxI[iz,:,:,:,2,1] * DxDiv +
      dXdxI[iz,:,:,:,2,2] * DyDiv) / J[iz,:,:,:]
  end  
end  

function GradDiv!(FuF,FvF,uF,vF,RhoF,Fe,dXdxI,J,Cache,Koeff)
  Nz = size(FuF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views Div = Cache.Block[:,:,:,1]
  @views temp = Cache.Block[:,:,:,2]
  @views DxDiv = Cache.Block[:,:,:,3]
  @views DyDiv = Cache.Block[:,:,:,4]

  @inbounds for iz = 1 : Nz
    @. Div = 0.0  
    @views @. temp = dXdxI[iz,:,:,:,1,1] * uF[iz,:,:,:] + dXdxI[iz,:,:,:,1,2] * vF[iz,:,:,:]
    DerivativeX!(Div,temp,DX)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * uF[iz,:,:,:] + dXdxI[iz,:,:,:,2,2] * vF[iz,:,:,:]
    DerivativeY!(Div,temp,DY)

    @. DxDiv = 0.0
    @. DyDiv = 0.0
    DerivativeX!(DxDiv,Div,DXW)
    DerivativeY!(DyDiv,Div,DYW)

    @views @. FuF[iz,:,:,:] -= Koeff * RhoF[iz,:,:,:] *(dXdxI[iz,:,:,:,1,1] * DxDiv +
      dXdxI[iz,:,:,:,2,1] * DyDiv) / J[iz,:,:,:]
    @views @. FvF[iz,:,:,:] -= Koeff * RhoF[iz,:,:,:] *(dXdxI[iz,:,:,:,2,1] * DxDiv +
      dXdxI[iz,:,:,:,2,2] * DyDiv) / J[iz,:,:,:]
  end  
end  

function RotCurl!(FuF,FvF,uF,vF,Fe,dXdxI,J,Cache)
  Nz = size(FuF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views W = Cache.Block[:,:,:,1]
  @views temp = Cache.Block[:,:,:,2]
  @views DxW = Cache.Block[:,:,:,2]
  @views DyW = Cache.Block[:,:,:,3]

  @inbounds for iz = 1 : Nz
    @. W = 0.0
    @views @. temp = dXdxI[iz,:,:,:,1,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,1,2] * uF[iz,:,:,:]
    DerivativeX!(W,temp,DX)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,2,2] * uF[iz,:,:,:]
    DerivativeY!(W,temp,DY)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DXW)
    DerivativeY!(DyW,W,DYW)
    @views @. FvF[iz,:,:,:] += (-dXdxI[iz,:,:,:,1,1] * DxW -
      dXdxI[iz,:,:,:,2,1] * DyW) / J[iz,:,:,:]
    @views @. FuF[iz,:,:,:] += (dXdxI[iz,:,:,:,2,1] * DxW +
      dXdxI[iz,:,:,:,2,2] * DyW) / J[iz,:,:,:]
  end   
end   

function RotCurl!(FuF,FvF,uF,vF,RhoF,Fe,dXdxI,J,Cache,Koeff)
  Nz = size(FuF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  DX = Fe.DX
  DY = Fe.DY
  DXW = Fe.DXW
  DYW = Fe.DYW

  @views W = Cache.Block[:,:,:,1]
  @views temp = Cache.Block[:,:,:,2]
  @views DxW = Cache.Block[:,:,:,2]
  @views DyW = Cache.Block[:,:,:,3]

  @inbounds for iz = 1 : Nz
    @. W = 0.0  
    @views @. temp = dXdxI[iz,:,:,:,1,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,1,2] * uF[iz,:,:,:]
    DerivativeX!(W,temp,DX)
    @views @. temp = dXdxI[iz,:,:,:,2,1] * vF[iz,:,:,:] - dXdxI[iz,:,:,:,2,2] * uF[iz,:,:,:]
    DerivativeY!(W,temp,DY)

    @. DxW = 0.0
    @. DyW = 0.0
    DerivativeX!(DxW,W,DXW)
    DerivativeY!(DyW,W,DYW)
    @views @. FvF[iz,:,:,:] = -Koeff * RhoF[iz,:,:,:] * (-dXdxI[iz,:,:,:,1,1] * DxW -
      dXdxI[iz,:,:,:,2,1] * DyW) / J[iz,:,:,:]
    @views @. FuF[iz,:,:,:] = -Koeff * RhoF[iz,:,:,:]  * (dXdxI[iz,:,:,:,2,1] * DxW +
      dXdxI[iz,:,:,:,2,2] * DyW) / J[iz,:,:,:]
  end
end
