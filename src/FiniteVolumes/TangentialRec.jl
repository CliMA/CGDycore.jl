mutable struct KiteFace
  MatTan::Array{Float64,2}
  LocGlob::Array{Int,1}
  KiteVol::Array{Float64,1}
end

function KiteFace()

  MatTan = zeros(0,0)
  LocGlob = zeros(Int,0)
  KiteVol = zeros(0)

  return KiteFace(
    MatTan,
    LocGlob,
    KiteVol,
  )
end  

mutable struct FaceT
 P1::Grids.Point
 P2::Grids.Point
 P3::Grids.Point
 P4::Grids.Point
end

function DivDiv(F,LocGlob,fDiv!,Ord,Rad,Jacobi)
  w,xw = GaussLobatto(Ord)
  W = zeros(length(w)^2,1)
  ksi = zeros(2,length(w)^2)
  iw = 1
  for i = 1 : length(w)
    for j = 1 : length(w)
      W[iw] = w[i] * w[j]
      ksi[1,iw] = xw[i]
      ksi[2,iw] = xw[j]
      iw=iw+1
    end
  end
  fDivLoc = zeros(4,length(W))
  for iw = 1 : length(W)
    @views fDiv!(fDivLoc[:,iw],ksi[1,iw],ksi[2,iw])
  end

  NumF = length(F)
  S=zeros(2*NumF,2*NumF)
  SLoc=zeros(4,4)
  for iF = 1 : NumF
    @. SLoc = 0
    for iw = 1 : length(W)
      dXdx,J,dXdxIT = Jacobi(F[iF],ksi[1,iw],ksi[2,iw],Rad)
      for i = 1 : 4
        for j = 1 : 4
          SLoc[i,j] = SLoc[i,j] + W[iw]*J*fDivLoc[i,iw] * fDivLoc[j,iw]
        end
      end
    end
    S[LocGlob[:,iF],LocGlob[:,iF]] = S[LocGlob[:,iF],LocGlob[:,iF]] + SLoc
  end
  return S
end

function RT0!(f,ksi1,ksi2)
  
  f[1,1] = 0
  f[2,1] = 1/2 * (1 - ksi2)

  f[1,2] = 1/2 * (1 + ksi1)
  f[2,2] = 0

  f[1,3] = 0
  f[2,3] = 1/2 * (1 + ksi2)

  f[1,4] = 1/2 * (1 - ksi1)
  f[2,4] = 0
end

function DivRT0!(f,ksi1,ksi2)
  f[1] = -1/2
  f[2] = 1/2
  f[3] = 1/2
  f[4] = -1/2
end

function Ned0!(f,ksi1,ksi2)

  f[1,1] = 1/2 * (1 - ksi2)
  f[2,1] = 0

  f[1,2] = 0
  f[2,2] = -1/2 * (1 + ksi1)

  f[1,3] = 1/2 * (1 + ksi2)
  f[2,3] = 0

  f[1,4] = 0
  f[2,4] = -1/2 * (1 - ksi1)
end

function MassMatrix(M,LocGlob,fCurl!,Ord)
  w,xw = GaussLobatto(Ord)
  W = zeros(length(w)^2,1)
  ksi = zeros(2,length(w)^2)
  iw = 1
  for i = 1 : length(w)
    for j = 1 : length(w)
      W[iw] = w[i] * w[j]
      ksi[1,iw] = xw[i]
      ksi[2,iw] = xw[j]
      iw=iw+1
    end
  end
  fCurlLoc = zeros(2,4,length(W))
  for iw = 1 : length(W)
    @views fCurl!(fCurlLoc[:,:,iw],ksi[1,iw],ksi[2,iw])
  end  
end

function CG1Grad!(f,ksi1,ksi2)
  #f(1) = 0.25*(1-ksi1)*(1-ksi2)

  f[1,1] =-0.25 * (1 - ksi2)
  f[2,1] =-0.25 * (1 - ksi1)
  #f(2)  = 0.25*(1+ksi1)*(1-ksi2)
  f[1,2] = 0.25 * (1 - ksi2)
  f[2,2] =-0.25 * (1 + ksi1)
  #f(3)  = 0.25*(1+ksi1)*(1+ksi2)
  f[1,3] = 0.25 * (1 + ksi2)
  f[2,3] = 0.25 * (1 + ksi1)
  #f(4)  = 0.25*(1-ksi1)*(1+ksi2)
  f[1,4] = -0.25 * (1 + ksi2)
  f[2,4] = 0.25 * (1 - ksi1)
  @. f = f * (3/4)
end

function WeakVort(F,LocGlob,IC,fGrad!,fu!,Ord,Rad,Jacobi)
  w,xw = GaussLobatto(Ord)
  W = zeros(length(w)^2,1)
  ksi = zeros(2,length(w)^2)
  iw = 1
  for i = 1 : length(w)
    for j = 1 : length(w)
      W[iw] = w[i] * w[j]
      ksi[1,iw] = xw[i]
      ksi[2,iw] = xw[j]
      iw=iw+1
    end
  end
  fGradLoc = zeros(2,4,length(W))
  fuLoc = zeros(2,4,length(W))
  for iw = 1 : length(W)
    @views fu!(fuLoc[:,:,iw],ksi[1,iw],ksi[2,iw])
    @views fGrad!(fGradLoc[:,:,iw],ksi[1,iw],ksi[2,iw])
  end
  NF = length(F)
  S  = zeros(1,2*NF)
  e = zeros(2*NF)
  for iF = 1 : NF
    e[1]=norm(F[iF].P1-F[iF].P2)
    e[2]=norm(F[iF].P2-F[iF].P3)
    e[3]=norm(F[iF].P3-F[iF].P4)
    e[4]=norm(F[iF].P4-F[iF].P1)

    SLoc = zeros(4)
    i = IC[iF]
    for iw = 1 : length(W)
      dXdx,J,dXdxIT,k = Jacobi(F[iF],ksi[1,iw],ksi[2,iw],Rad)
      uT = dXdxIT * fGradLoc[:,i,iw]
      # Rotate
      uTR = cross(uT,k)
      for j = 1 : 4
        SLoc[j] = SLoc[j] - W[iw] * uTR'*(dXdx*fuLoc[:,j,iw])*e[j]
      end
    end
    S[1,LocGlob[1:4,iF]] = S[1,LocGlob[1:4,iF]] + SLoc[:]
  end
  return S
end

function LocalGrid(Face,Grid)

  NumNodes = length(Face.N) 
  if NumNodes == 4
    F = Array{FaceT,1}(undef,4)
    P1 = Grid.Nodes[Face.N[1]].P
    P2 = Grid.Nodes[Face.N[2]].P
    P3 = Grid.Nodes[Face.N[3]].P
    P4 = Grid.Nodes[Face.N[4]].P
    PE1 = 0.5 * (P1 + P2)
    PE2 = 0.5 * (P2 + P3)
    PE3 = 0.5 * (P3 + P4)
    PE4 = 0.5 * (P4 + P1)
    PC = 0.25 * (P1 + P2 + P3 + P4)
    LocGlob = zeros(Int,4,4)
    IC = zeros(Int,4)

    F[1] = FaceT(P1,PE1,PC,PE4)
    LocGlob[1:4,1] = [1 5 8 4]
    IC[1]=3

    F[2] = FaceT(PE1,P2,PE2,PC)
    LocGlob[1:4,2] = [1 2 6 5]
    IC[2]=4

    F[3] = FaceT(PC,PE2,P3,PE3)
    LocGlob[1:4,3]=[6 2 3 7]
    IC[3]=1

    F[4] = FaceT(PE4,PC,PE3,P4)
    LocGlob[1:4,4]=[8 7 3 4]
    IC[4]=2
  elseif NumNodes == 5
    F = Array{FaceT,1}(undef,5)
    P1 = Grid.Nodes[Face.N[1]].P
    P2 = Grid.Nodes[Face.N[2]].P
    P3 = Grid.Nodes[Face.N[3]].P
    P4 = Grid.Nodes[Face.N[4]].P
    P5 = Grid.Nodes[Face.N[5]].P
    PE1 = 0.5 * (P1 + P2)
    PE2 = 0.5 * (P2 + P3)
    PE3 = 0.5 * (P3 + P4)
    PE4 = 0.5 * (P4 + P5)
    PE5 = 0.5 * (P5 + P1)
    PC = 1/5 * (P1 + P2 + P3 + P4 + P5)
    LocGlob = zeros(Int,4,5)
    IC = zeros(Int,5)

    F[1] = FaceT(P1,PE1,PC,PE5)
    LocGlob[1:4,1] = [1 6 10 5]
    IC[1]=3

    F[2] = FaceT(PE1,P2,PE2,PC)
    LocGlob[1:4,2] = [1 2 7 6]
    IC[2]=4

    F[3] = FaceT(PC,PE2,P3,PE3)
    LocGlob[1:4,3]=[7 2 3 8]
    IC[3]=1

    F[4] = FaceT(PC,PE3,P4,PE4)
    LocGlob[1:4,4]=[8 3 4 9]
    IC[4]=1

    F[5] = FaceT(PE5,PC,PE4,P5)
    LocGlob[1:4,5]=[10 9 4 5]
    IC[5]=2
  elseif NumNodes == 6
    F = Array{FaceT,1}(undef,6)
    P1 = Grid.Nodes[Face.N[1]].P
    P2 = Grid.Nodes[Face.N[2]].P
    P3 = Grid.Nodes[Face.N[3]].P
    P4 = Grid.Nodes[Face.N[4]].P
    P5 = Grid.Nodes[Face.N[5]].P
    P6 = Grid.Nodes[Face.N[6]].P
    PE1 = 0.5 * (P1 + P2)
    PE2 = 0.5 * (P2 + P3)
    PE3 = 0.5 * (P3 + P4)
    PE4 = 0.5 * (P4 + P5)
    PE5 = 0.5 * (P5 + P6)
    PE6 = 0.5 * (P6 + P1)
    PC = 1/6 * (P1 + P2 + P3 + P4 + P5 + P6)
    LocGlob = zeros(Int,4,6)
    IC = zeros(Int,6)

    F[1] = FaceT(P1,PE1,PC,PE6)
    LocGlob[1:4,1] = [1 7 12 6]
    IC[1] = 3

    F[2] = FaceT(PE1,P2,PE2,PC)
    LocGlob[1:4,2] = [1 2 8 7]
    IC[2] = 4

    F[3] = FaceT(PC,PE2,P3,PE3)
    LocGlob[1:4,3]=[8 2 3 9]
    IC[3] = 1

    F[4] = FaceT(PC,PE3,P4,PE4)
    LocGlob[1:4,4]=[9 3 4 10]
    IC[4] = 1

    F[5] = FaceT(PC,PE4,P5,PE5)
    LocGlob[1:4,5]=[10 5 5 11]
    IC[5] = 1

    F[6] = FaceT(PE6,PC,PE5,P6)
    LocGlob[1:4,6]=[12 11 5 6]
    IC[6] = 2

  end  
  return F,LocGlob,IC
end

function JacobiCart(F,ksi1,ksi2,Rad);
  dXdx = zeros(2,2)
  dXdx[1,1] = 0.25 * (-F.P1.x * (1 - ksi2) + F.P2.x * (1 - ksi2) + 
    F.P3.x * (1 + ksi2) - F.P4.x * (1+ksi2))
  dXdx[2,1] = 0.25 * (-F.P1.y * (1 - ksi2) + F.P2.y * (1 - ksi2) + 
    F.P3.y * (1 + ksi2) - F.P4.y * (1 + ksi2))

  dXdx[1,2] = 0.25 * (-F.P1.x * (1 - ksi1) - F.P2.x * (1 + ksi1) + 
    F.P3.x * (1 + ksi1) + F.P4.x * (1 - ksi1));
  dXdx[2,2] = 0.25 * (-F.P1.y * (1 - ksi1) - F.P2.y * (1 + ksi1) +
    F.P3.y * (1 + ksi1) + F.P4.y * (1 - ksi1));
  J = [dXdx 
       0.0 0.0]
  @views detJLoc = det(J[:,1],J[:,2])
  detJ = detJLoc
  pinvJ  = pinvJac(J)
  k = SVector{3}(0.0,0.0,1.0)
  return J,detJ,pinvJ,k
end


@inline function JacobiSphere(F,ksi1,ksi2,Rad)

  XT1 =  0.25*(F.P1.x*(1-ksi1)*(1-ksi2)+
            F.P2.x*(1+ksi1)*(1-ksi2)+
            F.P3.x*(1+ksi1)*(1+ksi2)+
            F.P4.x*(1-ksi1)*(1+ksi2))
    
  XT2 =  0.25*(F.P1.y*(1-ksi1)*(1-ksi2)+
            F.P2.y*(1+ksi1)*(1-ksi2)+
            F.P3.y*(1+ksi1)*(1+ksi2)+
            F.P4.y*(1-ksi1)*(1+ksi2))
           
  XT3 =  0.25*(F.P1.z*(1-ksi1)*(1-ksi2)+
           F.P2.z*(1+ksi1)*(1-ksi2)+
            F.P3.z*(1+ksi1)*(1+ksi2)+
            F.P4.z*(1-ksi1)*(1+ksi2))
    
  XLoc = SVector{3}(XT1,XT2,XT3)
  JP = @SArray[F.P1.x F.P2.x F.P3.x F.P4.x;
               F.P1.y F.P2.y F.P3.y F.P4.y;
               F.P1.z F.P2.z F.P3.z F.P4.z]

  J3 = @SArray([-0.25 + 0.25*ksi2  -0.25 + 0.25*ksi1
                 0.25 - 0.25*ksi2  -0.25 - 0.25*ksi1
                 0.25 + 0.25*ksi2   0.25 + 0.25*ksi1
                -0.25 - 0.25*ksi2   0.25 - 0.25*ksi1])


  f = Rad *(XT1^2 + XT2^2 + XT3^2)^(-3/2)
  dX1dXT1 = f * (XT2^2 + XT3^2)
  dX1dXT2= -f * XT1 * XT2 
  dX1dXT3= -f * XT1 * XT3
  dX2dXT1 = dX1dXT2 
  dX2dXT2 = f * (XT1^2+XT3^2)
  dX2dXT3 = -f * XT2 * XT3
  dX3dXT1 = dX1dXT3 
  dX3dXT2 = dX2dXT3 
  dX3dXT3 = f*(XT1^2+XT2^2)

  J1  =   @SArray([dX1dXT1    dX1dXT2     dX1dXT3   
                dX2dXT1     dX2dXT2     dX2dXT3
                dX3dXT1     dX3dXT2     dX3dXT3])   
  J   =   J1*JP*J3

  @views detJLoc = det(J[:,1],J[:,2])
  detJ = detJLoc
  pinvJ  = pinvJac(J)
  X = XLoc / norm(XLoc) 
  return J,detJ,pinvJ,X
end

@inline function det(a,b)

  d = (a[2] * b[3] - a[3] * b[2])^2 +
      (a[1] * b[3] - a[3] * b[1])^2 +
      (a[1] * b[2] - a[2] * b[1])^2
  d = sqrt(d)
end

@inline function pinvJac(J)
  g11 = J[1,1] * J[1,1] + J[2,1] * J[2,1] + J[3,1] * J[3,1]
  g12 = J[1,1] * J[1,2] + J[2,1] * J[2,2] + J[3,1] * J[3,2]
  g22 = J[1,2] * J[1,2] + J[2,2] * J[2,2] + J[3,2] * J[3,2]
  det = g11 * g22 - g12^2
  i11 = g22 / det
  i21 = -g12 / det
  i12 = -g12 / det
  i22 = g11 / det
  pJ11 = J[1,1] * i11 + J[1,2] * i21
  pJ12 = J[1,1] * i12 + J[1,2] * i22
  pJ21 = J[2,1] * i11 + J[2,2] * i21
  pJ22 = J[2,1] * i12 + J[2,2] * i22
  pJ31 = J[3,1] * i11 + J[3,2] * i21
  pJ32 = J[3,1] * i12 + J[3,2] * i22
  @SArray([pJ11 pJ12
          pJ21 pJ22
          pJ31 pJ32])

end

function GaussLobatto(Ord)
  w = zeros(Ord)
  xw = zeros(Ord)
  if Ord == 1
    w[1] = 2;
    xw[1] = 0
  elseif Ord == 2
    w[1] = 1
    w[2] = 1
    xw[1] = -1
    xw[2] = 1
  elseif Ord == 3
    w[1] = 1/3
    w[2] = 4/3
    w[3] = 1/3
    xw[1] = -1
    xw[2] = 0
    xw[3] = 1
  elseif Ord == 4
    w[1] = 1/6
    w[2] = 5/6
    w[3] = 5/6
    w[4] = 1/6
    xw[1] = -1
    xw[2] = -sqrt(1/5)
    xw[3] = sqrt(1/5)
    xw[4] = 1
  end
return w,xw
end

function HCurlHDiv()
end

function TangentialDiv(F,LocGlob,fu!,Ord,Rad,Jacobi)
  w,xw = GaussLobatto(Ord)
  fuLocB = zeros(2,4,length(w),4)
  for iw = 1 : length(w)
    @views fu!(fuLocB[:,:,iw,1],xw[iw],-1)
    @views fu!(fuLocB[:,:,iw,2],1,xw[iw])
    @views fu!(fuLocB[:,:,iw,3],xw[iw],1)
    @views fu!(fuLocB[:,:,iw,4],-1,xw[iw])
  end
  NumF = length(F)

  NVal = zeros(2*NumF,2*NumF)
  e = zeros(4)
  for iF = 1 : NumF
    SLoc=zeros(4,4)
    e[1] = Grids.norm(F[iF].P1-F[iF].P2)
    e[2] = Grids.norm(F[iF].P2-F[iF].P3)
    e[3] = Grids.norm(F[iF].P3-F[iF].P4)
    e[4] = Grids.norm(F[iF].P4-F[iF].P1)
    ss = 0.0
    for iw = 1 : length(w)
      # edge 1
      ie = 1
      t = [1;0]
      dXdx,J,dXdxIT = Jacobi(F[iF],xw[iw],-1,Rad)
      ss = ss + w[iw] * norm((dXdx * t))
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw]  * (dXdx * t)' * 
          (dXdx * fuLocB[:,j,iw,1]) / J
      end  
      # edge 2
      ie = 2
      t = [0;1]
      dXdx,J,dXdxIT = Jacobi(F[iF],1,xw[iw],Rad)
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw]  * (dXdx * t)' * 
          (dXdx * fuLocB[:,j,iw,2]) / J
      end
      # edge 3
      t = [1;0]
      ie = 3
      dXdx,J,dXdxIT = Jacobi(F[iF],xw[iw],1,Rad)
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw]  * (dXdx * t)' * 
          (dXdx * fuLocB[:,j,iw,3]) / J
      end
      t = [0;1]
      ie = 4
      dXdx,J,dXdxIT = Jacobi(F[iF],-1,xw[iw],Rad)
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw] * (dXdx * t)' *
          (dXdx * fuLocB[:,j,iw,4]) / J
       end
    end
    @show ss,e[1]
    NVal[LocGlob[:,iF],LocGlob[:,iF]] = NVal[LocGlob[:,iF],LocGlob[:,iF]] + SLoc
  end
  return NVal
end

function MatrixTangential(Jacobi,Grid)
  NumFaces = Grid.NumFaces

  KiteFaces=map(1:NumFaces) do i
    FiniteVolumes.KiteFace()
  end

  for iF = 1 : Grid.NumFaces
    F,LocGlob,IC = FiniteVolumes.LocalGrid(Grid.Faces[iF],Grid)
    NumF = length(F)

    Ord = 4;
    SDiv = FiniteVolumes.DivDiv(F,LocGlob,FiniteVolumes.DivRT0!,Ord,Grid.Rad,Jacobi);
    SDivV = FiniteVolumes.WeakVort(F,LocGlob,IC,FiniteVolumes.CG1Grad!,FiniteVolumes.RT0!,
      Ord,Grid.Rad,Jacobi);
    SSDiv=[SDiv;SDivV]

    for i = 1 : size(SSDiv,1)
      @show SSDiv[i,:]
    end  
    uB = zeros(NumF,NumF)
    for i = 1 : NumF
      uB[i,i] = .5;
    end
    uI = -SSDiv[NumF+1:end,NumF+1:end] \ (SSDiv[NumF+1:end,1:NumF] * uB)
    for i = 1 : size(uI,1)
      @show uI[i,:]
    end  
    NVal = TangentialDiv(F,LocGlob,FiniteVolumes.RT0!,Ord,Grid.Rad,Jacobi)
    for i = 1 : size(NVal,1)
      @show NVal[i,:]
    end  
    U = [uB;uI]
    T = NVal * U
    @. T[NumF+1:end,:] *= 0.5
    for i = 1 : size(T,1)
      @show T[i,:]
    end  
    KiteFaces[iF].MatTan = T
    LocGlob = zeros(Int,length(Grid.Faces[iF].E))
    @. LocGlob = Grid.Faces[iF].E
    KiteFaces[iF].LocGlob = LocGlob
    KiteVol = zeros(NumF)
    for i = 1 : NumF
      KiteVol[i] = (Grids.AreaSphericalTriangle(F[i].P1,F[i].P2,F[i].P3) + 
        Grids.AreaSphericalTriangle(F[i].P2,F[i].P3,F[i].P4)) * Grid.Rad^2
    end  
    KiteFaces[iF].KiteVol = KiteVol
    stop
  end    
  return KiteFaces
end    
