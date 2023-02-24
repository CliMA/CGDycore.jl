mutable struct FeElem
  OrdPolyX::Int
  OrdPolyZ::Int
  DX::Array{Float64, 2}
  DZT::Array{Float64, 2}
  IntZF2C::Array{Float64, 2}
  IntZC2F::Array{Float64, 2}
  wX::Array{Float64, 1}
  xw::Array{Float64, 1}
  xe::Array{Float64, 1}
  wZC::Array{Float64, 1}
  zwC::Array{Float64, 1}
  wZ::Array{Float64, 1}
  zw::Array{Float64, 1}
  ze::Array{Float64, 1}
  P::Array{Float64, 2}
  IntXF2cE::Array{Float64, 2}
  IntZC2cE::Array{Float64, 2}
  IntZF2cE::Array{Float64, 2}
  IntXE2F::Array{Float64, 2}
  IntZE2FT::Array{Float64, 2}
end

function FeElem(OrdPolyX,OrdPolyZ)
  #Horizontal grid
  (wX,xw)=GaussLobattoQuad(OrdPolyX)
  xe = -1.0 : 2.0/OrdPolyX : 1.0
  DX=zeros(OrdPolyX+1,OrdPolyX+1)
  for i=1:OrdPolyX+1
    for j=1:OrdPolyX+1
      DX[i,j]=DLagrange(xw[i],xw,j)
    end 
  end 
  #  Vertical Grid
  (wZ,zw)=GaussLobattoQuad(OrdPolyZ)
  (wZC,zwC)=GaussLegendreQuad(OrdPolyZ)
  ze = -1.0 : 2.0/OrdPolyZ : 1.0

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

  IntZE2FT = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for j = 1 : OrdPolyZ + 1
    for i = 1 : OrdPolyZ + 1
      IntZE2FT[j,i] = Lagrange(zw[i],ze,j) 
    end
  end  

  IntXF2cE = zeros(OrdPolyX,OrdPolyX+1)
  dx = 2.0 / OrdPolyX
  for j = 1 : OrdPolyX + 1
    xE = 0.5 * dx  
    for i = 1 : OrdPolyX
      IntXF2cE[i,j] = Lagrange(xE,xw,j) 
    xE = xE + dx
    end
  end  
  IntZC2cE = zeros(OrdPolyZ,OrdPolyZ)
  dz = 2.0 / OrdPolyZ
  for j = 1 : OrdPolyZ
    zE = 0.5 * dz  
    for i = 1 : OrdPolyZ 
      IntZC2cE[i,j] = Lagrange(zE,zwC,j) 
      zE = zE + dz
    end
  end  
  IntZF2cE = zeros(OrdPolyZ,OrdPolyZ+1)
  dz = 2.0 / OrdPolyZ
  for j = 1 : OrdPolyZ + 1
    zE = 0.5 * dz  
    for i = 1 : OrdPolyZ 
      IntZF2cE[i,j] = Lagrange(zE,zw,j) 
      zE = zE + dz
    end
  end  

  DZT = zeros(OrdPolyZ+1,OrdPolyZ+1)
  for i=1:OrdPolyZ+1
    for j=1:OrdPolyZ+1
      DZT[j,i]=DLagrange(zw[i],zw,j)
    end
  end
  P = (diagm(1.0./wZC) * IntZC2F') * diagm(wZ)

  return FeElem(
    OrdPolyX,
    OrdPolyZ,
    DX,
    DZT,
    IntZF2C,
    IntZC2F,
    wX,
    xw,
    xe,
    wZC,
    zwC,
    wZ,
    zw,
    ze,
    P,
    IntXF2cE,
    IntZC2cE,
    IntZF2cE,
    IntXE2F,
    IntZE2FT,
  )
end  
