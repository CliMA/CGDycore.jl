function TransSphere(ksi,eta,z,F,Topo,Param)
X=zeros(3)
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x)
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y)
X[3]=0.25*((1-ksi)*(1-eta)*F.P[1].z+
  (1+ksi)*(1-eta)*F.P[2].z+
  (1+ksi)*(1+eta)*F.P[3].z+
  (1-ksi)*(1+eta)*F.P[4].z)
r=norm(X)
X=X/r
X=X*(Param.RadPrint+z)
return X
end

function TransSphereS(ksi,eta,F,Rad)
X=zeros(3)
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x)
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y)
X[3]=0.25*((1-ksi)*(1-eta)*F.P[1].z+
  (1+ksi)*(1-eta)*F.P[2].z+
  (1+ksi)*(1+eta)*F.P[3].z+
  (1-ksi)*(1+eta)*F.P[4].z)
r=norm(X)
X=X/r
X=X*Rad
end

function TransSphereS!(P,ksi,eta,F,Rad)
  P.x=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x)
  P.y=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y)
  P.z=0.25*((1-ksi)*(1-eta)*F.P[1].z+
  (1+ksi)*(1-eta)*F.P[2].z+
  (1+ksi)*(1+eta)*F.P[3].z+
  (1-ksi)*(1+eta)*F.P[4].z)
  r=norm(P)
  P=P/r
  P=P*Rad
end

function TransSphereX(ksi,eta,zeta,X,CG,Global)
  OrdPoly=CG.OrdPoly
  OrdPolyZ=CG.OrdPolyZ
  XP=zeros(3)
  @inbounds for j=1:OrdPoly+1
    @inbounds for i=1:OrdPoly+1
      @inbounds for k=1:OrdPolyZ+1
        XP[:]=XP[:]+DG.Lagrange(ksi,CG.xw,i)*
          DG.Lagrange(eta,CG.xw,j)*DG.Lagrange(zeta,CG.xwZ,k)*X[i,j,k,:]
      end
    end
  end
  r=norm(XP)
  z=max(r-Global.Grid.Rad,0.0)
  XP=XP*(Global.Output.RadPrint+z)
  return XP
end

function TransSphereX!(XP,ksi,eta,zeta,X,CG,Global,::Grids.Quad)
  OrdPoly=CG.OrdPoly
  OrdPolyZ=CG.OrdPolyZ

  @. XP = 0
  iDoF = 1
  @inbounds for j = 1 : OrdPoly + 1 
    Lj = DG.Lagrange(eta,CG.xwCPU,j)
    @inbounds for i = 1 : OrdPoly + 1 
      Li = DG.Lagrange(ksi,CG.xwCPU,i) * Lj
      @inbounds for k = 1 : OrdPolyZ + 1 
        Fac = Li * DG.Lagrange(zeta,CG.xwZCPU,k)
        @inbounds for l = 1 : 3 
           XP[l] += Fac * X[iDoF,k,l]
        end 
      end 
      iDoF += 1
    end 
  end 

  r=norm(XP)
  z=max(r-Global.Grid.Rad,0.0)
  XP=XP/r
  XP=XP*(Global.Output.RadPrint+z)
end

function TransSphere2DX!(XP,ksi,eta,zeta,X,CG,Global,::Grids.Quad)
  OrdPoly=CG.OrdPoly

  @. XP = 0
  iDoF = 1
  @inbounds for j = 1 : OrdPoly + 1
    Lj = DG.Lagrange(eta,CG.xwCPU,j)
    @inbounds for i = 1 : OrdPoly + 1
      Fac = DG.Lagrange(ksi,CG.xwCPU,i) * Lj
      @inbounds for l = 1 : 3
         XP[l] += Fac * X[iDoF,l]
      end
      iDoF += 1
    end
  end

  r=norm(XP)
  z=max(r-Global.Grid.Rad,0.0)
  XP=XP/r
  XP=XP*(Global.Output.RadPrint+z)
end

function TransSphereX!(XP,ksi,eta,zeta,X,CG,Global,::Grids.Tri)
  DoF = CG.DoF
  OrdPolyZ=CG.OrdPolyZ
  @. XP = 0
  @inbounds for k = 1 : OrdPolyZ + 1
    XR = CG.PL2CPU * X[:,k,:]
    Fac = DG.Lagrange(zeta,CG.xwZCPU,k)
    @inbounds for i = 1 : size(XR,1)
      XP .+= Fac * CG.phi[i](ksi,eta) * XR[i,:]   
    end
  end  
  r=norm(XP)
  z=max(r-Global.Grid.Rad,0.0)
  XP=XP/r
  XP=XP*(Global.Output.RadPrint+z)
end

function TransCart(ksi,eta,z,F,Topo,Grid)
X=zeros(3);
X[1]=0.25*((1-ksi)*(1-eta)*F.P[1].x+
  (1+ksi)*(1-eta)*F.P[2].x+
  (1+ksi)*(1+eta)*F.P[3].x+
  (1-ksi)*(1+eta)*F.P[4].x);
X[2]=0.25*((1-ksi)*(1-eta)*F.P[1].y+
  (1+ksi)*(1-eta)*F.P[2].y+
  (1+ksi)*(1+eta)*F.P[3].y+
  (1-ksi)*(1+eta)*F.P[4].y);
X[3]=first(Topo(X[1],X[2],z,Grid.Topography));
return X
end

function TransCartX(ksi,eta,zeta,X,CG,Global)
  OrdPoly=CG.OrdPoly
  OrdPolyZ=CG.OrdPolyZ
  XP = zeros(3)
  @inbounds for j = 1 : OrdPoly + 1
    Lj = DG.Lagrange(eta,CG.xwCPU,j)
    @inbounds for i = 1 : OrdPoly + 1
      Li = DG.Lagrange(ksi,CG.xwCPU,i) * Lj
      @inbounds for k = 1 : OrdPolyZ + 1
        Fac = Li * DG.Lagrange(zeta,CG.xwZCPU,k)
        @inbounds for l = 1 : 3
           XP[l] += Fac * X[iDoF,k,l]
        end
      end
      iDoF += 1
    end
  end
  if abs(XP[1])<1.e-20
    XP[1] = 0.0  
  end  
  if abs(XP[2])<1.e-20
    XP[2] = 0.0  
  end  
  if abs(XP[3])<1.e-20
    XP[3] = 0.0  
  end  
  return XP
end  

function TransCartX!(XP,ksi,eta,zeta,X,CG,Global,::Grids.Quad)
  OrdPoly=CG.OrdPoly
  OrdPolyZ=CG.OrdPolyZ

  @. XP = 0
  iDoF = 1
  @inbounds for j = 1 : OrdPoly + 1
    Lj = DG.Lagrange(eta,CG.xwCPU,j)
    @inbounds for i = 1 : OrdPoly + 1
      Li = DG.Lagrange(ksi,CG.xwCPU,i) * Lj
      @inbounds for k = 1 : OrdPolyZ + 1
        Fac = Li * DG.Lagrange(zeta,CG.xwZCPU,k)
        @inbounds for l = 1 : 3
           XP[l] += Fac * X[iDoF,k,l]
        end
      end
      iDoF += 1
    end
  end

  if abs(XP[1])<1.e-20
    XP[1] = 0.0
  end
  if abs(XP[2])<1.e-20
    XP[2] = 0.0
  end
  if abs(XP[3])<1.e-20
    XP[3] = 0.0
  end
end

function TransCartX!(XP,ksi,eta,zeta,X,CG,Global,::Grids.Tri)
  DoF = CG.DoF
  OrdPolyZ=CG.OrdPolyZ
  @. XP = 0
  @inbounds for k = 1 : OrdPolyZ + 1
    XR = CG.PL2CPU * X[:,k,:]
    Fac = DG.Lagrange(zeta,CG.xwZCPU,k)
    @inbounds for i = 1 : size(XR,1)
      XP .+= Fac * CG.phi[i](ksi,eta) * XR[i,:]
    end
  end
  if abs(XP[1])<1.e-20
    XP[1] = 0.0
  end
  if abs(XP[2])<1.e-20
    XP[2] = 0.0
  end
  if abs(XP[3])<1.e-20
    XP[3] = 0.0
  end
end




