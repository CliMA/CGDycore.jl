function JacRiemannSparseH(DG,Metric,Grid,Aux,Phys)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]

  NH = Metric.NH
  VolSurfH = Metric.VolSurfH
  wF = DG.wF
  NE = size(DG.GlobE,3)
  nz = size(Metric.dXdxI,5)
  NF = size(Metric.dXdxI,6)
  N = size(DG.DVT,1)
  M = size(DG.DVZT,1)

  @views Theta = Aux[:,:,:,4]
  @views pRhoTh = Aux[:,:,:,3]
  @show NE,N
  for iE = 1 : NE
    iFL = Grid.EF[1,iE]
    iFR = Grid.EF[2,iE]
    for i = 1 : N
      indL = DG.GlobE[1,i,iE]
      indR = DG.GlobE[2,i,iE]
      for iz = 1 : nz
        for k = 1 : M  
          n1 = NH[k,iz,i,iE,1]
          n2 = NH[k,iz,i,iE,2]
          n3 = NH[k,iz,i,iE,3]  
          Surf = VolSurfH[k,iz,i,iE] / wF[i]

          # Rho
          iRowL = Index(k,indL,iz,RhoPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,RhoPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * Surf
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLoc,RowInd,ColInd,Val)
          # v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * Surf
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLoc,RowInd,ColInd,Val)
          # w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * Surf
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLoc,RowInd,ColInd,Val)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * Phys.invcS * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * Phys.invcS * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLocR,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLocL,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLocR,RowInd,ColInd,Val)

          # RhoTh
          iRowL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.25 * n1 * Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLoc,RowInd,ColInd,Val)
          # v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.25 * n2 * Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLoc,RowInd,ColInd,Val)
          # w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.25 * n3 * Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLoc,RowInd,ColInd,Val)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.25 * Phys.invcS * pRhoTh[k,iz,indL] * 
            Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          ValLocR = 0.25 * Phys.invcS * pRhoTh[k,iz,indR] * 
            Surf * (Theta[k,iz,indL] + Theta[k,iz,indR])
          Insert!(iRowL,iColL,-ValLocL,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLocR,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLocL,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLocR,RowInd,ColInd,Val)

          # u
          iRowL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,uPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * n1 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          #v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * n2 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          #w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n1 * n3 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * n1 * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * n1 * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLocR,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLocL,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLocR,RowInd,ColInd,Val)

          # v
          iRowL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,vPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * n1 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          #v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * n2 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          #w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n2 * n3 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * n2 * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * n2 * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLocR,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLocL,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLocR,RowInd,ColInd,Val)

          # w
          iRowL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iRowR = Index(k,indR,iz,wPos,M,nz,N,NF)
          # u
          iColL = Index(k,indL,iz,uPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,uPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * n1 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          #v
          iColL = Index(k,indL,iz,vPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,vPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * n2 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          #w
          iColL = Index(k,indL,iz,wPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,wPos,M,nz,N,NF)
          ValLoc = 0.5 * n3 * n3 * Surf * Phys.cS
          Insert!(iRowL,iColL,-ValLoc,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLoc,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,-ValLoc,RowInd,ColInd,Val)
          # RhoTh
          iColL = Index(k,indL,iz,RhoThPos,M,nz,N,NF)
          iColR = Index(k,indR,iz,RhoThPos,M,nz,N,NF)
          ValLocL = 0.5 * n3 * pRhoTh[k,iz,indL] * Surf
          ValLocR = 0.5 * n3 * pRhoTh[k,iz,indR] * Surf
          Insert!(iRowL,iColL,-ValLocL,RowInd,ColInd,Val)
          Insert!(iRowL,iColR,-ValLocR,RowInd,ColInd,Val)
          Insert!(iRowR,iColL,ValLocL,RowInd,ColInd,Val)
          Insert!(iRowR,iColR,ValLocR,RowInd,ColInd,Val)
        end
      end  
    end
  end  
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  return Jac
end

function JacRiemannSparseV(DG,Metric,Aux,Phys)

  RhoPos = 1
  uPos = 2
  vPos = 3
  wPos = 4
  RhoThPos = 5

  RowInd = Int[]
  ColInd = Int[]
  Val = Float64[]

  NVG = Metric.NV
  VolSurfVG = Metric.VolSurfV
  nz = size(Metric.dXdxI,5)
  N = size(DG.DVT,1)
  M = size(DG.DVZT,1)
  NF = size(Metric.dXdxI,6)
  VolSurfV = reshape(VolSurfVG,nz+1,N,N,NF)
  NV = reshape(NVG,nz+1,N,N,NF,3)
  Theta = reshape(Aux[:,:,:,4],M,nz,N,N,NF)
  pRhoTh = reshape(Aux[:,:,:,3],M,nz,N,N,NF)

  for iF = 1 : NF
    for j = 1 : N
      for i = 1 : N
        for iz = 1 : nz + 1
          n1 = NV[iz,i,j,iF,1]
          n2 = NV[iz,i,j,iF,2]
          n3 = NV[iz,i,j,iF,3]  
          Surf = VolSurfV[iz,i,j,iF]
          if iz == 1
            # u
            iRowT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            # u
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n1 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # v
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n2 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # w
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n3 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n1 * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)

            # v
            iRowT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            # u
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n1 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # v
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n2 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # w
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n3 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n2 * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)

            # w
            iRowT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            # u
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n1 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # v
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n2 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # w
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n3 * Surf
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n3 * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
          elseif iz == nz + 1
            # u
            iRowB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n1 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n1 * pRhoTh[M,iz-1,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)

            # v
            iRowB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n2 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n2 * pRhoTh[M,iz-1,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)

            # w
            iRowB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            ValLoc = Phys.cS * n3 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            ValLoc = n3 * pRhoTh[M,iz-1,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
          else
            # Rho 
            iRowB = Index(M,i,j,iz-1,iF,RhoPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,RhoPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * Surf
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * Phys.invcS * pRhoTh[M,iz-1,i,j,iF] * Surf
            ValLocT = 0.5 * Phys.invcS * pRhoTh[1,iz,i,j,iF] * Surf
            Insert!(iRowB,iColB,-ValLocB,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLocT,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLocB,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLocT,RowInd,ColInd,Val)


            # RhoTh 
            iRowB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.25 * n1 * Surf * (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF])
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.25 * n2 * Surf * (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF])
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.25 * n3 * Surf * (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF])
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.25 * Phys.invcS * pRhoTh[M,iz-1,i,j,iF] * 
              (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF]) * Surf
            ValLocT = 0.25 * Phys.invcS * pRhoTh[1,iz,i,j,iF] * 
              (Theta[M,iz-1,i,j,iF] + Theta[1,iz,i,j,iF]) * Surf
            Insert!(iRowB,iColB,-ValLocB,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLocT,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLocB,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLocT,RowInd,ColInd,Val)

            # u 
            iRowB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * n1 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * n2 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n1 * n3 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * n1 * pRhoTh[M,iz-1,i,j,iF] * Surf 
            ValLocT = 0.5 * n1 * pRhoTh[1,iz,i,j,iF] * Surf 
            Insert!(iRowB,iColB,-ValLocB,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLocT,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLocB,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLocT,RowInd,ColInd,Val)

            # v 
            iRowB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * n1 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * n2 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n2 * n3 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * n2 * pRhoTh[M,iz-1,i,j,iF] * Surf 
            ValLocT = 0.5 * n2 * pRhoTh[1,iz,i,j,iF] * Surf 
            Insert!(iRowB,iColB,-ValLocB,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLocT,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLocB,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLocT,RowInd,ColInd,Val)

            # w 
            iRowB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iRowT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            # u
            iColB = Index(M,i,j,iz-1,iF,uPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,uPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * n1 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # v
            iColB = Index(M,i,j,iz-1,iF,vPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,vPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * n2 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # w
            iColB = Index(M,i,j,iz-1,iF,wPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,wPos,M,nz,N,NF)
            ValLoc = 0.5 * n3 * n3 * Surf * Phys.cS
            Insert!(iRowB,iColB,-ValLoc,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLoc,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,-ValLoc,RowInd,ColInd,Val)
            # RhoTh
            iColB = Index(M,i,j,iz-1,iF,RhoThPos,M,nz,N,NF)
            iColT = Index(1,i,j,iz,iF,RhoThPos,M,nz,N,NF)
            ValLocB = 0.5 * n3 * pRhoTh[M,iz-1,i,j,iF] * Surf 
            ValLocT = 0.5 * n3 * pRhoTh[1,iz,i,j,iF] * Surf 
            Insert!(iRowB,iColB,-ValLocB,RowInd,ColInd,Val)
            Insert!(iRowB,iColT,-ValLocT,RowInd,ColInd,Val)
            Insert!(iRowT,iColB,ValLocB,RowInd,ColInd,Val)
            Insert!(iRowT,iColT,ValLocT,RowInd,ColInd,Val)
          end
        end
      end
    end
  end  
  NZe = M * nz * N * N * NF * 5
  Jac = sparse(RowInd, ColInd, Val, NZe, NZe)
  return Jac
end

