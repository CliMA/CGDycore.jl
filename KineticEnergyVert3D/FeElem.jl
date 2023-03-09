mutable struct FeElem
  OrdPolyX::Int
  OrdPolyY::Int
  OrdPolyZ::Int
  DX::Array{Float64, 2}
  DXW::Array{Float64, 2}
  DY::Array{Float64, 2}
  DYW::Array{Float64, 2}
  DZ::Array{Float64, 2}
  DZT::Array{Float64, 2}
  IntZF2C::Array{Float64, 2}
  IntZC2F::Array{Float64, 2}
  wX::Array{Float64, 1}
  xw::Array{Float64, 1}
  xe::Array{Float64, 1}
  wY::Array{Float64, 1}
  yw::Array{Float64, 1}
  ye::Array{Float64, 1}
  wZC::Array{Float64, 1}
  zwC::Array{Float64, 1}
  wZ::Array{Float64, 1}
  zw::Array{Float64, 1}
  ze::Array{Float64, 1}
  P::Array{Float64, 2}
  IntXF2cE::Array{Float64, 2}
  IntYF2cE::Array{Float64, 2}
  IntZC2cE::Array{Float64, 2}
  IntZF2cE::Array{Float64, 2}
  IntXE2F::Array{Float64, 2}
  IntYE2F::Array{Float64, 2}
  IntZE2F::Array{Float64, 2}
  IntZE2FT::Array{Float64, 2}
end

function FeElem(OrdPolyX,OrdPolyY,OrdPolyZ)
  #Horizontal grid
  (wX,xw)=GaussLobattoQuad(OrdPolyX)
  xe = zeros(OrdPolyX+1)
  xe[1] = -1.0
  for i = 2 : OrdPolyX
    xe[i] = xe[i-1] + 2.0/OrdPolyX
  end  
  xe[OrdPolyX+1] = 1.0
  DX=zeros(OrdPolyX+1,OrdPolyX+1)
  for i=1:OrdPolyX+1
    for j=1:OrdPolyX+1
      DX[i,j]=DLagrange(xw[i],xw,j)
    end 
  end 
  DXW = -inv(diagm(wX))*DX'*diagm(wX)

  (wY,yw)=GaussLobattoQuad(OrdPolyY)
  ye = zeros(OrdPolyY+1)
  ye[1] = -1.0
  for i = 2 : OrdPolyY
    ye[i] = ye[i-1] + 2.0/OrdPolyY
  end  
  ye[OrdPolyY+1] = 1.0
  DY=zeros(OrdPolyY+1,OrdPolyY+1)
  for i=1:OrdPolyY+1
    for j=1:OrdPolyY+1
      DY[i,j]=DLagrange(yw[i],yw,j)
    end 
  end 
  DYW = -inv(diagm(wY))*DY'*diagm(wY)
  #  Vertical Grid
  (wZ,zw)=GaussLobattoQuad(OrdPolyZ)
  (wZC,zwC)=GaussLegendreQuad(OrdPolyZ)
  ze = zeros(OrdPolyZ+1)
  ze[1] = -1.0
  for i = 2 : OrdPolyZ
    ze[i] = ze[i-1] + 2.0/OrdPolyZ
  end  
  ze[OrdPolyZ+1] = 1.0

  IntZC2F = zeros(Float64,OrdPolyZ+1,OrdPolyZ)
  IntZF2C = zeros(Float64,OrdPolyZ,OrdPolyZ+1)
  for j = 1 : OrdPolyZ+1
    for i = 1 : OrdPolyZ
      IntZF2C[i,j]=Lagrange(zwC[i],zw,j)
      IntZC2F[j,i]=Lagrange(zw[j],zwC,i)
    end 
  end 

  IntXE2F = zeros(OrdPolyX+1,OrdPolyX+1)
  for j = 1 : OrdPolyX + 1
    for i = 1 : OrdPolyX +1
      IntXE2F[i,j] = Lagrange(xw[i],xe,j) 
    end
  end  
  IntYE2F = zeros(OrdPolyY+1,OrdPolyY+1)
  for j = 1 : OrdPolyY + 1
    for i = 1 : OrdPolyY +1
      IntYE2F[i,j] = Lagrange(yw[i],ye,j) 
    end
  end  

  IntZE2FT = zeros(OrdPolyZ+1,OrdPolyZ+1)
  IntZE2F = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for j = 1 : OrdPolyZ + 1
    for i = 1 : OrdPolyZ + 1
      IntZE2F[i,j] = Lagrange(zw[i],ze,j) 
      IntZE2FT[j,i] = Lagrange(zw[i],ze,j) 
    end
  end  

  IntXF2cE = zeros(OrdPolyX,OrdPolyX+1)
  dx = 2.0 / OrdPolyX
  for j = 1 : OrdPolyX + 1
    xE = -1.0 + 0.5 * dx  
    for i = 1 : OrdPolyX
      IntXF2cE[i,j] = Lagrange(xE,xw,j) 
      xE = xE + dx
    end
  end  
  IntYF2cE = zeros(OrdPolyY,OrdPolyY+1)
  dy = 2.0 / OrdPolyY
  for j = 1 : OrdPolyY + 1
    yE = -1.0 + 0.5 * dy  
    for i = 1 : OrdPolyY
      IntYF2cE[i,j] = Lagrange(yE,yw,j) 
    yE = yE + dy
    end
  end  
  IntZC2cE = zeros(OrdPolyZ,OrdPolyZ)
  dz = 2.0 / OrdPolyZ
  for j = 1 : OrdPolyZ
    zE = -1.0 + 0.5 * dz  
    for i = 1 : OrdPolyZ 
      IntZC2cE[i,j] = Lagrange(zE,zwC,j) 
      zE = zE + dz
    end
  end  
  IntZF2cE = zeros(OrdPolyZ,OrdPolyZ+1)
  dz = 2.0 / OrdPolyZ
  for j = 1 : OrdPolyZ + 1
    zE = -1.0 + 0.5 * dz  
    for i = 1 : OrdPolyZ 
      IntZF2cE[i,j] = Lagrange(zE,zw,j) 
      zE = zE + dz
    end
  end  

  DZ = zeros(OrdPolyZ+1,OrdPolyZ+1)
  DZT = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for i=1:OrdPolyZ+1
    for j=1:OrdPolyZ+1
      DZ[i,j]=DLagrange(zw[i],zw,j)
      DZT[j,i]=DLagrange(zw[i],zw,j)
    end
  end
  P = (diagm(1.0./wZC) * IntZC2F') * diagm(wZ)

  return FeElem(
    OrdPolyX,
    OrdPolyY,
    OrdPolyZ,
    DX,
    DXW,
    DY,
    DYW,
    DZ,
    DZT,
    IntZF2C,
    IntZC2F,
    wX,
    xw,
    xe,
    wY,
    yw,
    ye,
    wZC,
    zwC,
    wZ,
    zw,
    ze,
    P,
    IntXF2cE,
    IntYF2cE,
    IntZC2cE,
    IntZF2cE,
    IntXE2F,
    IntYE2F,
    IntZE2F,
    IntZE2FT,
  )
end  

function DerivativeX!(Dc,c,DX)
  nx = size(c,1)
  ny = size(c,2)
  nz = size(c,3)
  @inbounds for k = 1 : nz
    @inbounds for j = 1 : ny
      @inbounds for i = 1 : nx
        @inbounds for l = 1 : nx
          Dc[i,j,k] = Dc[i,j,k] + DX[i,l] * c[l,j,k]
        end
      end
    end
  end
end

function DerivativeY!(Dc,c,DY)
  nx = size(c,1)
  ny = size(c,2)
  nz = size(c,3)
  @inbounds for k = 1 : nz
    @inbounds for i = 1 : nx
      @inbounds for j = 1 : ny
        @inbounds for l = 1 : ny
          Dc[i,j,k] = Dc[i,j,k] + DY[j,l] * c[i,l,k]
        end
      end
    end
  end
end

function DerivativeZ!(Dc,c,DZ)
  nx = size(c,1)
  ny = size(c,2)
  nz = size(c,3)
  @inbounds for j = 1 : ny
    @inbounds for i = 1 : nx
      @inbounds for k = 1 : nz
        @inbounds for l = 1 : nz
          Dc[i,j,k] = Dc[i,j,k] + DZ[k,l] * c[i,j,l]
        end
      end
    end
  end
end
