mutable struct KiteFace
  MatTan::Array{Float64,2}
  LocGlob::Array{Int,1}
end

function KiteFace()

  MatTan = zeros(0,0)
  LocGlob = zeros(Int,0)

  return KiteFace(
    MatTan,
    LocGlob,
  )
end  

mutable struct FaceT
 P1::Grids.Point
 P2::Grids.Point
 P3::Grids.Point
 P4::Grids.Point
end

function DivDiv(F,LocGlob,fDiv!,Ord)
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
      dXdx,J,dXdxIT = Jacobi(F[iF],ksi[1,iw],ksi[2,iw])
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

function WeakVort(F,LocGlob,IC,fGrad!,fu!,Ord)
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
      dXdx,J,dXdxIT = Jacobi(F[iF],ksi[1,iw],ksi[2,iw])
      uT = dXdxIT * fGradLoc[:,i,iw]
      # Rotate
      uTR = [uT[2];-uT[1]]
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
    LocGlob = zeros(Int,4,4)
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

function Jacobi(F,ksi1,ksi2);
  dXdx = zeros(2,2)
  dXdx[1,1] = 0.25 * (-F.P1.x * (1 - ksi2) + F.P2.x * (1 - ksi2) + 
    F.P3.x * (1 + ksi2) - F.P4.x * (1+ksi2))
  dXdx[2,1] = 0.25 * (-F.P1.y * (1 - ksi2) + F.P2.y * (1 - ksi2) + 
    F.P3.y * (1 + ksi2) - F.P4.y * (1 + ksi2))

  dXdx[1,2] = 0.25 * (-F.P1.x * (1 - ksi1) - F.P2.x * (1 + ksi1) + 
    F.P3.x * (1 + ksi1) + F.P4.x * (1 - ksi1));
  dXdx[2,2] = 0.25 * (-F.P1.y * (1 - ksi1) - F.P2.y * (1 + ksi1) +
    F.P3.y * (1 + ksi1) + F.P4.y * (1 - ksi1));
  J = det(dXdx)
  dXdxIT = inv(dXdx')
  return dXdx,J,dXdxIT
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

function TangentialDiv(F,LocGlob,fu!,Ord)
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
    for iw = 1 : length(w)
      # edge 1
      ie = 1
      t = [1;0]
      dXdx,J,dXdxIT = Jacobi(F[iF],xw[iw],-1)
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw] / e[1] * (dXdx * t)' * 
          (dXdx * fuLocB[:,j,iw,1]) * e[j]
      end  
      # edge 2
      ie = 2
      t = [0;1]
      dXdx,J,dXdxIT = Jacobi(F[iF],1,xw[iw])
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw] / e[2] * (dXdx * t)' * 
          (dXdx * fuLocB[:,j,iw,2]) * e[j]
      end
      # edge 3
      t = [1;0]
      ie = 3
      dXdx,J,dXdxIT = Jacobi(F[iF],xw[iw],1)
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw] / e[3] * (dXdx * t)' * 
          (dXdx * fuLocB[:,j,iw,3]) * e[j]
      end
      t = [0;1]
      ie = 4
      dXdx,J,dXdxIT = Jacobi(F[iF],-1,xw[iw])
      for j = 1 : 4
        SLoc[ie,j] = SLoc[ie,j] + w[iw] / e[4] * (dXdx * t)' *
          (dXdx * fuLocB[:,j,iw,4]) * e[j]
       end
    end
    NVal[LocGlob[:,iF],LocGlob[:,iF]] = NVal[LocGlob[:,iF],LocGlob[:,iF]] + SLoc
  end
  return NVal
end

function MatrixTangential(Grid)
  NumFaces = Grid.NumFaces

  KiteFaces=map(1:NumFaces) do i
    FiniteVolumes.KiteFace()
  end

  for iF = 1 : Grid.NumFaces
    F,LocGlob,IC = FiniteVolumes.LocalGrid(Grid.Faces[1],Grid)
    NumF = length(F)

    Ord = 4;
    SDiv = FiniteVolumes.DivDiv(F,LocGlob,FiniteVolumes.DivRT0!,Ord);
    SDivV = FiniteVolumes.WeakVort(F,LocGlob,IC,FiniteVolumes.CG1Grad!,FiniteVolumes.RT0!,Ord);
    SSDiv=[SDiv;SDivV]

    uB = zeros(NumF,NumF)
    for i = 1 : NumF
      uB[i,i] = 1;
    end
    uI = -SSDiv[NumF+1:end,NumF+1:end] \ (SSDiv[NumF+1:end,1:NumF] * uB)
    @show uI
    NVal = FiniteVolumes.TangentialDiv(F,LocGlob,FiniteVolumes.RT0!,Ord)
    @show NVal
    U = [uB;uI]
    T = NVal * U
    KiteFaces[iF].MatTan = T
    LocGlob = zeros(Int,length(Grid.Faces[iF].E))
    @. LocGlob = Grid.Faces[iF].E
    KiteFaces[iF].LocGlob = LocGlob
    stop
  end    
  return KiteFaces
end    
