function MetricVolume3(backend,FT,Grid)
  Faces = Grid.Faces
  Edges = Grid.Edges
  Nodes = Grid.Nodes
  NumFaces = Grid.NumFaces
  NumEdges = Grid.NumEdges
  NumNodes = Grid.NumNodes
  nz = Grid.nz
  Rad = Grid.Rad

  PrimalVolume = KernelAbstractions.zeros(backend,FT,nz,NumFaces)
  PrimalMidPoints = KernelAbstractions.zeros(backend,FT,3,nz,NumFaces)
  PrimalSideFaces = KernelAbstractions.zeros(backend,FT,nz,NumEdges)
  PrimalSideNormals = KernelAbstractions.zeros(backend,FT,3,nz,NumEdges)
  PrimalSideMidPoints = KernelAbstractions.zeros(backend,FT,3,nz,NumEdges)
  PrimalTopFaces = KernelAbstractions.zeros(backend,FT,nz+1,NumFaces)
  PrimalTopNormals = KernelAbstractions.zeros(backend,FT,3,nz+1,NumFaces)
  PrimalTopMidPoints = KernelAbstractions.zeros(backend,FT,3,nz+1,NumFaces)
  PrimalPoints = KernelAbstractions.zeros(backend,FT,3,nz+1,NumNodes)

  NQH = 5
  (xwH,wH)=gausslobatto(NQH)
  NQV = 5
  (xwV,wV)=gausslobatto(NQV)
  DH = zeros(NQH,NQH)
  for i = 1 : NQH
    for j = 1 : NQH
      DH[i,j] = DG.DLagrange(xwH[i],xwH,j)
      if abs(DH[i,j]) <= 1.e-12
        DH[i,j] = 0.0
      end
    end
  end
  DV = zeros(NQV,NQV)
  for i = 1 : NQV
    for j = 1 : NQV
      DV[i,j] = DG.DLagrange(xwV[i],xwV,j)
      if abs(DV[i,j]) <= 1.e-12
        DV[i,j] = 0.0
      end
    end
  end

  FiniteElements.EdgeFace!(backend,Grid)
  F = FiniteElements.PointsFromGrid(backend,FT,Grid)
# FillX!(backend,Metric,FE,F,Grid,zS)
# FillContravariant!(backend,Metric,FE,Grid,Grid.Type,Model.MetricType)
# FillDet!(backend,Metric,FE,Grid)
  X = zeros(NQH,NQH,NQV,3)
  dXdxI = zeros(3,3)
  AdaptGrid = Grids.AdaptGrid(FT,"Sleve",Grid.H)
  BottomNormal = zeros(3)
  TopNormal = zeros(3)
  SideXLNormal = zeros(3)
  SideXRNormal = zeros(3)
  SideYLNormal = zeros(3)
  SideYRNormal = zeros(3)
  CellMP2 = zeros(3)
  CellMP3 = zeros(3)
  CellMP = zeros(3)
  BottomMP = zeros(3)
  TopMP = zeros(3)
  SideXLMP = zeros(3)
  SideXRMP = zeros(3)
  SideYLMP = zeros(3)
  SideYRMP = zeros(3)
  for iF = 1 : Grid.NumFaces
    for iz = 1 : Grid.nz
      z1 = Grid.z[iz]
      z2 = Grid.z[iz+1]
      for j = 1 : NQH
        for i = 1 : NQH  
          for k = 1 : NQV  
            @views FiniteElements.XPoint!(AdaptGrid,X[i,j,k,:],Rad,
              xwH[i],xwH[j],xwV[k],F[:,:,iF],z1,z2,Grid.H,0.0,Grid.Type,Grid.Form)
          end  
        end  
      end  
      Vol = 0.0
      @. CellMP = 0.0
      AreaBottom = 0.0
      AreaTop = 0.0
      AreaSideXL = 0.0
      AreaSideXR = 0.0
      AreaSideYL = 0.0
      AreaSideYR = 0.0
      @. BottomNormal = 0.0
      @. TopNormal = 0.0
      @. SideXLNormal = 0.0
      @. SideXRNormal = 0.0
      @. SideYLNormal = 0.0
      @. SideYRNormal = 0.0
      @. BottomMP = 0.0
      @. TopMP = 0.0
      @. SideXLMP = 0.0
      @. SideXRMP = 0.0
      @. SideYLMP = 0.0
      @. SideYRMP = 0.0
      for j = 1 : NQH
        Vol2 = 0.0  
        @. CellMP2 = 0.0
        for i = 1 : NQH  
          Vol3 = 0.0  
          @. CellMP3 = 0.0
          for k = 1 : NQV  
            DXdx = DH[i,1] * X[1,j,k,1]
            DYdx = DH[i,1] * X[1,j,k,2]
            DZdx = DH[i,1] * X[1,j,k,3]
            DXdy = DH[j,1] * X[i,1,k,1]
            DYdy = DH[j,1] * X[i,1,k,2]
            DZdy = DH[j,1] * X[i,1,k,3]
            for l = 2 : NQH
              DXdx += DH[i,l] * X[l,j,k,1]
              DYdx += DH[i,l] * X[l,j,k,2]
              DZdx += DH[i,l] * X[l,j,k,3]
              DXdy += DH[j,l] * X[i,l,k,1]
              DYdy += DH[j,l] * X[i,l,k,2]
              DZdy += DH[j,l] * X[i,l,k,3]
            end
            DXdz = DV[k,1] * X[i,j,1,1]
            DYdz = DV[k,1] * X[i,j,1,2]
            DZdz = DV[k,1] * X[i,j,1,3]
            for l = 2 : NQV
              DXdz += DV[k,l] * X[i,j,l,1]
              DYdz += DV[k,l] * X[i,j,l,2]
              DZdz += DV[k,l] * X[i,j,l,3]
            end
            dXdxI[1,1] = DYdy * DZdz - DZdy * DYdz
            dXdxI[1,2] = DZdy * DXdz - DXdy * DZdz
            dXdxI[1,3] = DXdy * DYdz - DYdy * DXdz

            dXdxI[2,1] = DYdz * DZdx - DZdz * DYdx
            dXdxI[2,2] = DZdz * DXdx - DXdz * DZdx
            dXdxI[2,3] = DXdz * DYdx - DYdz * DXdx

            dXdxI[3,1] = DYdx * DZdy - DZdx * DYdy
            dXdxI[3,2] = DZdx * DXdy - DXdx * DZdy
            dXdxI[3,3] = DXdx * DYdy - DYdx * DXdy
            J = sqrt(dXdxI[1,1] * (dXdxI[2,2] * dXdxI[3,3] - dXdxI[2,3] * dXdxI[3,2]) -
              dXdxI[1,2] * (dXdxI[2,1] * dXdxI[3,3] - dXdxI[2,3] * dXdxI[3,1]) +
              dXdxI[1,3] * (dXdxI[2,1] * dXdxI[3,2] - dXdxI[2,2] * dXdxI[3,1]))
            Vol3 += J * wV[k]
            @. @views CellMP3  += X[i,j,k,:] * J * wV[k]
            if k == 1
              temp = sqrt(dXdxI[3,1]^2 + dXdxI[3,2]^2 + dXdxI[3,3]^2) *  wH[i] * wH[j]
              AreaBottom += temp
              @. BottomNormal += dXdxI[3,:] * wH[i] * wH[j]
              @. BottomMP += X[i,j,k,:] * temp
            end    
            if k == NQV
              temp = sqrt(dXdxI[3,1]^2 + dXdxI[3,2]^2 + dXdxI[3,3]^2) *  wH[i] * wH[j]
              AreaTop += temp
              @. @views TopNormal += dXdxI[3,:] * wH[i] * wH[j]
              @. @views TopMP += X[i,j,k,:] * temp
            end    
            if i == 1
              temp = sqrt(dXdxI[1,1]^2 + dXdxI[1,2]^2 + dXdxI[1,3]^2) *  wV[k] * wH[j]  
              AreaSideXL += temp
              @. @views SideXLNormal += dXdxI[1,:] *  wV[k] * wH[j]
              @. @views SideXLMP += X[i,j,k,:] * temp
            end  
            if i == NQH
              temp = sqrt(dXdxI[1,1]^2 + dXdxI[1,2]^2 + dXdxI[1,3]^2) *  wV[k] * wH[j]  
              AreaSideXR += temp
              @. @views SideXRNormal += dXdxI[1,:] * wV[k] * wH[j]
              @. @views SideXRMP += X[i,j,k,:] * temp
            end  
            if j == 1
              temp = sqrt(dXdxI[2,1]^2 + dXdxI[2,2]^2 + dXdxI[2,3]^2) *  wV[k] * wH[i]  
              AreaSideYL += temp
              @. @views SideYLNormal += dXdxI[2,:] * wV[k] * wH[i]
              @. @views SideYLMP += X[i,j,k,:] * temp
            end  
            if j == NQH
              temp = sqrt(dXdxI[2,1]^2 + dXdxI[2,2]^2 + dXdxI[2,3]^2) *  wV[k] * wH[i]  
              AreaSideYR += temp
              @. @views SideYRNormal += dXdxI[2,:] * wV[k] * wH[i]
              @. @views SideYRMP += X[i,j,k,:] * temp
            end  
          end  
          Vol2 += Vol3 * wH[i]
          @. CellMP2  += CellMP3 * wH[i]
        end  
        Vol += Vol2 * wH[j]
        @. CellMP  += CellMP2 * wH[j]
      end  
      PrimalVolume[iz,iF] = Vol
      NormNormal = norm(BottomNormal)
      @. @views PrimalMidPoints[:,iz,iF] = CellMP / Vol
      PrimalTopFaces[iz,iF] = AreaBottom
      @. @views  PrimalTopNormals[:,iz,iF] = BottomNormal / NormNormal
      @. @views  PrimalTopMidPoints[:,iz,iF] = BottomMP / AreaBottom

      NormNormal = norm(TopNormal)
      PrimalTopFaces[iz+1,iF] = AreaTop
      @. @views PrimalTopNormals[:,iz+1,iF] = TopNormal / NormNormal
      @. @views  PrimalTopMidPoints[:,iz+1,iF] = TopMP / AreaTop

      iE1 = Grid.Faces[iF].E[1]
      NormNormal = norm(SideYLNormal)
      PrimalSideFaces[iz,iE1] = AreaSideYL
      @. @views PrimalSideNormals[:,iz,iE1] = SideYLNormal / NormNormal
      @. @views PrimalSideMidPoints[:,iz,iE1] = SideYLMP / AreaSideYL

      iE2 = Grid.Faces[iF].E[2]
      NormNormal = norm(SideXRNormal)
      PrimalSideFaces[iz,iE2] = AreaSideXR
      @. @views PrimalSideNormals[:,iz,iE2] = SideXRNormal / NormNormal
      @. @views PrimalSideMidPoints[:,iz,iE2] = SideXRMP / AreaSideXR

      iE3 = Grid.Faces[iF].E[3]
      NormNormal = norm(SideYRNormal)
      PrimalSideFaces[iz,iE3] = AreaSideYR
      @. @views PrimalSideNormals[:,iz,iE3] = SideYRNormal / NormNormal
      @. @views PrimalSideMidPoints[:,iz,iE3] = SideYRMP / AreaSideYR

      iE4 = Grid.Faces[iF].E[4]
      NormNormal = norm(SideXLNormal)
      PrimalSideFaces[iz,iE4] = AreaSideXL
      @. @views PrimalSideNormals[:,iz,iE4] = SideXLNormal / NormNormal
      @. @views PrimalSideMidPoints[:,iz,iE4] = SideXLMP / AreaSideXL

      iN1 = Grid.Faces[iF].N[1]
      @. @views PrimalPoints[:,iz,iN1] = X[1,1,1,:]
      @. @views PrimalPoints[:,iz+1,iN1] = X[1,1,NQV,:]
      iN2 = Grid.Faces[iF].N[2]
      @. @views PrimalPoints[:,iz,iN2] = X[NQH,1,1,:]
      @. @views PrimalPoints[:,iz+1,iN2] = X[NQH,1,NQV,:]
      iN3 = Grid.Faces[iF].N[3]
      @. @views PrimalPoints[:,iz,iN3] = X[NQH,NQH,1,:]
      @. @views PrimalPoints[:,iz+1,iN3] = X[NQH,NQH,NQV,:]
      iN4 = Grid.Faces[iF].N[4]
      @. @views PrimalPoints[:,iz,iN4] = X[1,NQH,1,:]
      @. @views PrimalPoints[:,iz+1,iN4] = X[1,NQH,NQV,:]
    end  
  end  
  return MetricFiniteVolume3{FT,
                    typeof(PrimalVolume), 
                    typeof(PrimalSideNormals)}(
    PrimalVolume,
    PrimalMidPoints,
    PrimalSideFaces,
    PrimalTopFaces,
    PrimalSideNormals,
    PrimalTopNormals,
    PrimalSideMidPoints,
    PrimalTopMidPoints,
    PrimalPoints,
    )
end  




