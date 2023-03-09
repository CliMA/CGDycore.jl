function Advection!(FC,c,uC,vC,Fe,Metric,Cache)

  Nx = size(FC,1)
  Ny = size(FC,2)
  Nx = size(FC,3)

  @views @.  FC = 0.0
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        @views AdvectionBlock!(FC[ix,iy,iz,:,:],uC[ix,iy,iz,:,:],vC[ix,iy,iz,:,:],Fe,
          Metric.dXdxI,Metric.J,Cache)
      end
    end
  end  
  @views @. JRhoC = (Metric.J[:,:,:,:,:,1] + Metric.J[:,:,:,:,:,2]) * RhoC
  DSSC!(FC,JRhoC)
end

function AdvectionBlock!(Fc,uC,vC,Fe,dXdxI,J,Cache)
  DXW = Fe.DXW
  DYW = Fe.DYW
  @views W = Cache.Block[:,:,1,1]
  @views temp = Cache.Block[:,:,1,2]
  @views DxW = Cache.Block[:,:,1,2]
  @views DyW = Cache.Block[:,:,1,3]

  @views @. temp = 0.5 * c * ((dXdxI[:,:,1,1,1] + dXdxI[:,:,2,1,1]) * uC +
    (dXdxI[:,:,1,1,2] + dXdxI[:,:,2,1,2]) * vC)
  DerivtiveX!(Fc,temp,DXW)
  @views @. temp = 0.5 * c * ((dXdxI[:,:,1,2,1] + dXdxI[:,:,2,2,1]) * uC +
    (dXdxI[:,:,1,2,2] + dXdxI[:,:,2,2,2]) * vC)
  DerivtiveY!(Fc,temp,DYW)
end  

