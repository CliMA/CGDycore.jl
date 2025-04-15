mutable struct vtkStruct{FT<:AbstractFloat,
                         AT6<:AbstractArray}
  vtkInter::AT6
  cells
  pts::Array{Float64,2}
  RefineMidPoints::Array{Float64,2}
end
function vtkStruct{FT}(backend) where FT<:AbstractFloat
  vtkInter = KernelAbstractions.zeros(backend,FT,0,0,0,0,0,0)
  pts = Array{Float64,2}(undef,0,0)
  cells = MeshCell[]
  RefineMidPoints = zeros(0,0)
  return vtkStruct{FT,
                   typeof(vtkInter)}(
  vtkInter,
  cells,
  pts,
  RefineMidPoints,
  )
end

function vtkStruct{FT}(backend,Grid,NumFaces,Flat;Refine=0) where FT<:AbstractFloat

  vtkInter = KernelAbstractions.zeros(backend,FT,0,0,0,0,0,0)
  cells = MeshCell[]
  dTol = 2*pi / 30
  NumNodes = 0

  if Refine == 0
    celltype = VTKCellTypes.VTK_POLYGON
    for iF in 1 : NumFaces
      inds = Vector(1 : length(Grid.Faces[iF].N)) .+ NumNodes
      push!(cells, MeshCell(celltype, inds))
      NumNodes += length(Grid.Faces[iF].N)
    end
    pts = Array{Float64,2}(undef,3,NumNodes)
    NumNodes = 0
    RefineMidPoints = zeros(1,2)
    if Grid.Type == Grids.Tri()
      RefineMidPoints[1,1] = -1/3  
      RefineMidPoints[1,2] = -1/3  
    end  
    for iF in 1 : NumFaces
      if Grid.Form == "Sphere"
        if Flat
          lam = zeros(length(Grid.Faces[iF].N))
          theta = zeros(length(Grid.Faces[iF].N))
          NumNodesLoc = 0
          for iN in Grid.Faces[iF].N
            NumNodesLoc += 1
            (lam[NumNodesLoc],theta[NumNodesLoc],z) = Grids.cart2sphere(Grid.Nodes[iN].P.x,
              Grid.Nodes[iN].P.y,Grid.Nodes[iN].P.z)
          end
          lammin = minimum(lam)
          lammax = maximum(lam)
          if abs(lammin - lammax) > 2*pi-dTol
            for i = 1 : NumNodesLoc
              if lam[i] > pi
                lam[i] = lam[i] - 2*pi
                if lam[i] > 3*pi
                  lam[i] = lam[i]  - 2*pi
                end
              end
            end
          end
          for i = 1 : NumNodesLoc
            NumNodes += 1
            pts[:,NumNodes] = [lam[i],theta[i],0.0]
          end
        else    
          for iN in  Grid.Faces[iF].N
            NumNodes += 1
            pts[1,NumNodes] = Grid.Nodes[iN].P.x
            pts[2,NumNodes] = Grid.Nodes[iN].P.y
            pts[3,NumNodes] = Grid.Nodes[iN].P.z
          end
        end  
      else
        for i = 1 : length(Grid.Faces[iF].N)
          NumNodes += 1
          pts[1,NumNodes] = Grid.Faces[iF].P[i].x
          pts[2,NumNodes] = Grid.Faces[iF].P[i].y
          pts[3,NumNodes] = Grid.Faces[iF].P[i].z
        end
      end
    end
  else
    if Grid.Type == Grids.Quad()
      NodeLoc = zeros(3,4)
      lam = zeros(4)
      theta = zeros(4)
      celltypeQ = VTKCellTypes.VTK_QUAD
      NumRefine = (Refine + 1) * (Refine + 1)  
      RefineMidPoints = zeros(NumRefine,2)
      RefinePoints = zeros(NumRefine,4,2)
      dksi = 2 / (Refine + 1)
      iRefine = 0
      for jR in 1 : Refine + 1
        for iR in 1 : Refine + 1
          iRefine += 1  
          RefinePoints[iRefine,1,1] = -1.0 + (iR - 1) * dksi  
          RefinePoints[iRefine,1,2] = -1.0 + (jR - 1) * dksi  
          RefinePoints[iRefine,2,1] = -1.0 + (iR    ) * dksi  
          RefinePoints[iRefine,2,2] = -1.0 + (jR - 1) * dksi  
          RefinePoints[iRefine,3,1] = -1.0 + (iR    ) * dksi  
          RefinePoints[iRefine,3,2] = -1.0 + (jR    ) * dksi  
          RefinePoints[iRefine,4,1] = -1.0 + (iR - 1) * dksi  
          RefinePoints[iRefine,4,2] = -1.0 + (jR    ) * dksi  
          RefineMidPoints[iRefine,1] = 0.25 * (RefinePoints[iRefine,1,1] +
            RefinePoints[iRefine,2,1] + RefinePoints[iRefine,3,1] +
            RefinePoints[iRefine,4,1])
          RefineMidPoints[iRefine,2] = 0.25 * (RefinePoints[iRefine,1,2] +
            RefinePoints[iRefine,2,2] + RefinePoints[iRefine,3,2] +
            RefinePoints[iRefine,4,2])
        end
      end  
      for iF in 1 : NumFaces
        for i in 1 : NumRefine
          inds = Vector(1 : 4) .+ NumNodes
          push!(cells, MeshCell(celltypeQ, inds))
          NumNodes += 4
        end
      end  
      pts = Array{Float64,2}(undef,3,NumNodes)
      NumNodes = 0
      for iF in 1 : NumFaces
        for iR in 1 : NumRefine  
          if Grid.Form == "Sphere" 
            for i in 1 : 4  
              NodeLoc[1,i], NodeLoc[2,i], NodeLoc[3,i] =
                Bilinear(RefinePoints[iR,i,1],RefinePoints[iR,i,2],
                Grid.Nodes[Grid.Faces[iF].N[1]].P,
                Grid.Nodes[Grid.Faces[iF].N[2]].P,
                Grid.Nodes[Grid.Faces[iF].N[3]].P,
                Grid.Nodes[Grid.Faces[iF].N[4]].P,
                Grid.Form,Grid.Rad)
            end  
            if Flat
              for i in 1 : 3
                (lam[i],theta[i],z) = Grids.cart2sphere(NodeLoc[1,i],NodeLoc[2,i],NodeLoc[3,i])
              end
              lammin = minimum(lam)
              lammax = maximum(lam)
              if abs(lammin - lammax) > 2*pi-dTol
                for i = 1 : 4
                  if lam[i] > pi
                    lam[i] = lam[i] - 2*pi
                    if lam[i] > 3*pi
                      lam[i] = lam[i]  - 2*pi
                    end
                  end
                end
              end
              for i in 1 : 4
                NumNodes += 1
                pts[:,NumNodes] = [lam[i],theta[i],0.0]
              end
            else
              for i in 1 : 4
                NumNodes += 1
                pts[1,NumNodes] = NodeLoc[1,i]
                pts[2,NumNodes] = NodeLoc[2,i]
                pts[3,NumNodes] = NodeLoc[3,i]
              end
            end
          else    
            for i in 1 : 4  
              NumNodes += 1  
              pts[1,NumNodes], pts[2,NumNodes], pts[3,NumNodes] =
                Bilinear(RefinePoints[iR,i,1],RefinePoints[iR,i,2],
                Grid.Faces[iF].P[1],Grid.Faces[iF].P[2],Grid.Faces[iF].P[3],
                Grid.Faces[iF].P[4],Grid.Form,Grid.Rad)
            end    
          end    
        end    
      end    
    elseif Grid.Type == Grids.Tri()
      NodeLoc = zeros(3,3)
      lam = zeros(3)
      theta = zeros(3)
      celltypeQ = VTKCellTypes.VTK_TRIANGLE
      NumRefine = (Refine + 1) * (Refine + 1)  
      RefineMidPoints = zeros(NumRefine,2)
      RefinePoints = zeros(NumRefine,4,2)
      dksi = 2 / (Refine + 1)
      iRefine = 0
      for jRR in Refine + 1 : -1 : 1
        jR = Refine + 2 - jRR  
        for iR in 1 : jRR
          iRefine += 1  
          RefinePoints[iRefine,1,1] = -1.0 + (iR - 1) * dksi  
          RefinePoints[iRefine,1,2] = -1.0 + (jR - 1) * dksi  
          RefinePoints[iRefine,2,1] = -1.0 + (iR    ) * dksi  
          RefinePoints[iRefine,2,2] = -1.0 + (jR - 1) * dksi  
          RefinePoints[iRefine,3,1] = -1.0 + (iR - 1) * dksi  
          RefinePoints[iRefine,3,2] = -1.0 + (jR    ) * dksi  
          RefineMidPoints[iRefine,1] = 1/3 * (RefinePoints[iRefine,1,1] +
            RefinePoints[iRefine,2,1] + RefinePoints[iRefine,3,1]) 
          RefineMidPoints[iRefine,2] = 1/3 * (RefinePoints[iRefine,1,2] +
            RefinePoints[iRefine,2,2] + RefinePoints[iRefine,3,2])
          if iR == jRR
            continue
          end  
          iRefine += 1  
          RefinePoints[iRefine,1,1] = -1.0 + (iR    ) * dksi  
          RefinePoints[iRefine,1,2] = -1.0 + (jR - 1) * dksi  
          RefinePoints[iRefine,2,1] = -1.0 + (iR    ) * dksi  
          RefinePoints[iRefine,2,2] = -1.0 + (jR    ) * dksi  
          RefinePoints[iRefine,3,1] = -1.0 + (iR - 1) * dksi  
          RefinePoints[iRefine,3,2] = -1.0 + (jR    ) * dksi  
          RefineMidPoints[iRefine,1] = 1/3 * (RefinePoints[iRefine,1,1] +
            RefinePoints[iRefine,2,1] + RefinePoints[iRefine,3,1]) 
          RefineMidPoints[iRefine,2] = 1/3 * (RefinePoints[iRefine,1,2] +
            RefinePoints[iRefine,2,2] + RefinePoints[iRefine,3,2])
        end  
      end  
      for iF in 1 : NumFaces
        for i in 1 : NumRefine
          inds = Vector(1 : 3) .+ NumNodes
          push!(cells, MeshCell(celltypeQ, inds))
          NumNodes += 3
        end
      end  
      pts = Array{Float64,2}(undef,3,NumNodes)
      NumNodes = 0
      for iF in 1 : NumFaces
        for iR in 1 : NumRefine  
          if Grid.Form == "Sphere" 
            for i in 1 : 3  
              NodeLoc[1,i], NodeLoc[2,i], NodeLoc[3,i] =
                Linear(RefinePoints[iR,i,1],RefinePoints[iR,i,2],
                Grid.Nodes[Grid.Faces[iF].N[1]].P,
                Grid.Nodes[Grid.Faces[iF].N[2]].P,
                Grid.Nodes[Grid.Faces[iF].N[3]].P,
                Grid.Form,Grid.Rad)
            end    
            if Flat 
              for i in 1 : 3
                (lam[i],theta[i],z) = Grids.cart2sphere(NodeLoc[1,i],NodeLoc[2,i],NodeLoc[3,i])
              end  
              lammin = minimum(lam)
              lammax = maximum(lam)
              if abs(lammin - lammax) > 2*pi-dTol
                for i = 1 : 3
                  if lam[i] > pi
                    lam[i] = lam[i] - 2*pi
                    if lam[i] > 3*pi
                      lam[i] = lam[i]  - 2*pi
                    end
                  end
                end
              end
              for i in 1 : 3
                NumNodes += 1
                pts[:,NumNodes] = [lam[i],theta[i],0.0]
              end
            else    
              for i in 1 : 3  
                NumNodes += 1  
                pts[1,NumNodes] = NodeLoc[1,i]
                pts[2,NumNodes] = NodeLoc[2,i]
                pts[3,NumNodes] = NodeLoc[3,i]
              end
            end  
          else    
            for i in 1 : 3  
              NumNodes += 1  
              pts[1,NumNodes], pts[2,NumNodes], pts[3,NumNodes] =
                Linear(RefinePoints[iR,i,1],RefinePoints[iR,i,2],
                Grid.Faces[iF].P[1],Grid.Faces[iF].P[2],Grid.Faces[iF].P[3],
                Grid.Form,Grid.Rad)
            end    
          end    
        end    
      end    
    end
  end
  return vtkStruct{FT,
                   typeof(vtkInter)}(
  vtkInter,
  cells,
  pts,
  RefineMidPoints,
  )
end                   

#=
function vtkStruct{FT}(backend,OrdPrint::Int,Trans,FE::DyCore.FEStruct,Metric,Global) where FT<:AbstractFloat
  OrdPoly = FE.OrdPoly
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  Npts = 8 * NF * nz * OrdPrint * OrdPrint
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  x = zeros(8,3)
  lam=zeros(8,1)
  theta=zeros(8,1)
  z=zeros(8,1)
  dTol = Global.Output.dTol
  
  FTX = eltype(Metric.X)
  X = zeros(FTX,OrdPoly+1,OrdPoly+1,2,3)
  for iF = 1 : NF
    for iz = 1 : nz
      @views copyto!(X,reshape(Metric.X[:,:,:,iz,iF],OrdPoly+1,OrdPoly+1,2,3))
      dd = 2 / OrdPrint
      eta0 = -1
      for jRef = 1 : OrdPrint
        ksi0 = -1
        eta1 = eta0 + dd
        for iRef = 1 : OrdPrint
          ksi1 = ksi0 + dd
          @views Trans(x[1,:],ksi0,eta0, -1.0,X,FE,Global)
          @views Trans(x[2,:],ksi1,eta0, -1.0,X,FE,Global)
          @views Trans(x[3,:],ksi1,eta1, -1.0,X,FE,Global)
          @views Trans(x[4,:],ksi0,eta1, -1.0,X,FE,Global)
          @views Trans(x[5,:],ksi0,eta0, 1.0,X,FE,Global)
          @views Trans(x[6,:],ksi1,eta0, 1.0,X,FE,Global)
          @views Trans(x[7,:],ksi1,eta1, 1.0,X,FE,Global)
          @views Trans(x[8,:],ksi0,eta1, 1.0,X,FE,Global)
          if Global.Grid.Form == "Sphere" && Global.Output.Flat
            for i=1:8
              (lam[i],theta[i],z[i]) = Grids.cart2sphere(x[i,1],x[i,2],x[i,3])
            end 
            lammin = minimum(lam)
            lammax = maximum(lam)
            if abs(lammin - lammax) > 2*pi-dTol
              for i = 1 : 8
                if lam[i] > pi
                  lam[i] = lam[i] - 2*pi
                  if lam[i] > 3*pi
                    lam[i] = lam[i]  - 2*pi
                  end
                end
              end
            end
            for i = 1 : 8
              pts[:,ipts] = [lam[i],theta[i],max(z[i]-Global.Grid.Rad,0.0)/Global.Grid.H*3]
              ipts = ipts + 1
            end
          else
            for i=1:8
              pts[:,ipts] = [x[i,1],x[i,2],x[i,3]]
              ipts = ipts + 1
            end
          end
          ksi0=ksi1
        end
        eta0=eta1
      end
    end
  end
  celltype = VTKCellTypes.VTK_HEXAHEDRON

  ConnectivityList=reshape(1:1:8*NF*OrdPrint*OrdPrint*nz,8,NF*OrdPrint*OrdPrint*nz)
  cells = MeshCell[]

  for k in 1 : NF * nz * OrdPrint * OrdPrint
    inds = Vector(1 : 8) .+ 8 * (k -1)
    push!(cells, MeshCell(celltype, inds))
  end

  vtkInter = zeros(Float64,OrdPrint,OrdPrint,1,OrdPoly+1,OrdPoly+1,1)
  dd=2/OrdPrint
  eta0=-1
  for jRef=1:OrdPrint
    ksi0=-1
    eta1=eta0+dd
    for iRef=1:OrdPrint
      ksi1=ksi0+dd
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          vtkInter[iRef,jRef,1,i,j,1] = vtkInter[iRef,jRef,i,j] + DG.Lagrange(0.5*(ksi0+ksi1),FE.xwCPU,i)*
                DG.Lagrange(0.5*(eta0+eta1),FE.xwCPU,j)
        end
      end
      ksi0 = ksi1
    end
    eta0 = eta1
  end
  dvtkInter = KernelAbstractions.zeros(backend,FT,size(vtkInter))
  copyto!(dvtkInter,vtkInter)
  return vtkStruct{FT,
                   typeof(dvtkInter)}(
    dvtkInter,
    cells,
    pts,
  )  
end
=#

function vtkStruct{FT}(backend,OrdPrint::Int,Trans,FE,Metric,Global) where FT<:AbstractFloat
  OrdPoly = FE.OrdPoly
  OrdPolyZ = FE.OrdPolyZ
  OrdPrintZ = Global.Output.OrdPrintZ
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  Npts = 8 * NF * nz * OrdPrint * OrdPrint * OrdPrintZ
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  x = zeros(8,3)
  lam=zeros(8,1)
  theta=zeros(8,1)
  z=zeros(8,1)
  dTol = Global.Output.dTol
  dTol=2*pi/max(Global.Output.nPanel-1,2)
  
  FTX = eltype(Metric.X)
  X = zeros(FTX,OrdPoly+1,OrdPoly+1,OrdPolyZ+1,3)
  dd = 2 / OrdPrint
  ddZ = 2 / OrdPrintZ
  for iF = 1 : NF
    for iz = 1 : nz
      @views copyto!(X,reshape(Metric.X[:,:,:,iz,iF],OrdPoly+1,OrdPoly+1,OrdPolyZ+1,3))
      zeta0 = -1.0
      for kRef = 1 : OrdPrintZ  
        zeta1 = zeta0 + ddZ
        eta0 = -1.0
        for jRef = 1 : OrdPrint
          ksi0 = -1.0
          eta1 = eta0 + dd
          for iRef = 1 : OrdPrint
            ksi1 = ksi0 + dd
            @views Trans(x[1,:],ksi0,eta0,zeta0,X,FE,Global)
            @views Trans(x[2,:],ksi1,eta0,zeta0,X,FE,Global)
            @views Trans(x[3,:],ksi1,eta1,zeta0,X,FE,Global)
            @views Trans(x[4,:],ksi0,eta1,zeta0,X,FE,Global)
            @views Trans(x[5,:],ksi0,eta0,zeta1,X,FE,Global)
            @views Trans(x[6,:],ksi1,eta0,zeta1,X,FE,Global)
            @views Trans(x[7,:],ksi1,eta1,zeta1,X,FE,Global)
            @views Trans(x[8,:],ksi0,eta1,zeta1,X,FE,Global)
            if Global.Grid.Form == "Sphere" && Global.Output.Flat
              for i=1:8
                (lam[i],theta[i],z[i]) = Grids.cart2sphere(x[i,1],x[i,2],x[i,3])
              end 
              lammin = minimum(lam)
              lammax = maximum(lam)
              if abs(lammin-lammax) > pi
                smaller = count(<(pi), lam)
                greater = count(>=(pi), lam)

                if greater >= smaller
                  for i = 1 : 8
                    if lam[i] < pi
                      lam[i] += 2 * pi
                    end
                  end
                else
                  for i = 1 : 8
                    if lam[i] >= pi
                      lam[i] -= 2 * pi
                    end
                  end
                end
              end
          #=
              if abs(lammin - lammax) > 2*pi-dTol
                @show iF,lam[1:4]  
                for i = 1 : 8
                  if lam[i] > 2*pi-dTol
                    lam[i] = lam[i] - 2*pi
                  end  
                  if lam[i] < dTol
                    lam[i] = lam[i]  + 2*pi
                  end
                end
                @show iF,lam[1:4]  
              end
          =#    
              for i = 1 : 8
                pts[:,ipts] = [lam[i],theta[i],max(z[i]-Global.Grid.Rad,0.0)/Global.Grid.H*3]
                ipts = ipts + 1
              end
            else
              for i=1:8
                pts[:,ipts] = [x[i,1],x[i,2],x[i,3]]
                ipts = ipts + 1
              end
            end
            ksi0=ksi1
          end
          eta0=eta1
        end
        zeta0=zeta1
      end
    end
  end
  celltype = VTKCellTypes.VTK_HEXAHEDRON

  ConnectivityList=reshape(1:1:8*NF*OrdPrint*OrdPrint*OrdPrintZ*nz,8,NF*OrdPrint*OrdPrint*OrdPrintZ*nz)
  cells = MeshCell[]

  for k in 1 : NF * nz * OrdPrint * OrdPrint * OrdPrintZ
    inds = Vector(1 : 8) .+ 8 * (k -1)
    push!(cells, MeshCell(celltype, inds))
  end

  if typeof(FE) <: FiniteElements.CGQuad
    vtkInter = zeros(Float64,OrdPrint,OrdPrint,1,OrdPoly+1,OrdPoly+1,1)
    eta0=-1.0
    zeta1=zeta0+ddZ
    for jRef=1:OrdPrint
      ksi0=-1.0
      eta1=eta0+dd
      for iRef=1:OrdPrint
        ksi1=ksi0+dd
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            vtkInter[iRef,jRef,1,i,j,1] = DG.Lagrange(0.5*(ksi0+ksi1),FE.xwCPU,i)*
                DG.Lagrange(0.5*(eta0+eta1),FE.xwCPU,j) 
          end
        end
        ksi0 = ksi1
      end
      eta0 = eta1
    end
  else    
    vtkInter = zeros(Float64,OrdPrint,OrdPrint,OrdPrintZ,OrdPoly+1,OrdPoly+1,OrdPolyZ+1)
    zeta0 = -1.0
    for kRef = 1 : OrdPrintZ
      eta0=-1.0
      zeta1=zeta0+ddZ
      for jRef=1:OrdPrint
        ksi0=-1.0
        eta1=eta0+dd
        for iRef=1:OrdPrint
          ksi1=ksi0+dd
          for k=1:OrdPolyZ+1
            for j=1:OrdPoly+1
              for i=1:OrdPoly+1
                vtkInter[iRef,jRef,kRef,i,j,k] = DG.Lagrange(0.5*(ksi0+ksi1),FE.xwCPU,i)*
                    DG.Lagrange(0.5*(eta0+eta1),FE.xwCPU,j) *
                    DG.Lagrange(0.5*(zeta0+zeta1),FE.xwZCPU,k) 
              end
            end
          end
          ksi0 = ksi1
        end
        eta0 = eta1
      end
      zeta0 = zeta1
    end
  end
  dvtkInter = KernelAbstractions.zeros(backend,FT,size(vtkInter))
  copyto!(dvtkInter,vtkInter)
  RefineMidPoints = zeros(0,0)
  return vtkStruct{FT,
                   typeof(dvtkInter)}(
    dvtkInter,
    cells,
    pts,
    RefineMidPoints,
  )  
end

function vtkInit2D(OrdPrint::Int,Trans,FE,Metric,Global)
  OrdPoly = FE.OrdPoly
  NF = Global.Grid.NumFaces
  Npts = 4 * NF * OrdPrint * OrdPrint
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  x = zeros(4,3)
  lam=zeros(4,1)
  theta=zeros(4,1)
  z=zeros(4,1)
  dTol = Global.Output.dTol
  FT = eltype(Metric.X)
  backend = get_backend(Metric.X)
  X = zeros(FT,OrdPoly+1,OrdPoly+1,1,3)
  for iF = 1 : NF
    @views copyto!(X,reshape(Metric.X[:,1,:,1,iF],OrdPoly+1,OrdPoly+1,1,3))
    dd = 2 / OrdPrint
    eta0 = -1
    for jRef = 1 : OrdPrint
      ksi0 = -1
      eta1 = eta0 + dd
      for iRef = 1 : OrdPrint
        ksi1 = ksi0 + dd
        @views Trans(x[1,:],ksi0,eta0, -1.0,X,FE,Global)
        @views Trans(x[2,:],ksi1,eta0, -1.0,X,FE,Global)
        @views Trans(x[3,:],ksi1,eta1, -1.0,X,FE,Global)
        @views Trans(x[4,:],ksi0,eta1, -1.0,X,FE,Global)
        if Global.Grid.Form == "Sphere" && Global.Output.Flat
          for i=1:4
            (lam[i],theta[i],z[i]) = Grids.cart2sphere(x[i,1],x[i,2],x[i,3])
          end 
          lammin = minimum(lam)
          lammax = maximum(lam)
          if abs(lammin - lammax) > 2*pi-dTol
            for i = 1 : 4
              if lam[i] > pi
                lam[i] = lam[i] - 2*pi
                if lam[i] > 3*pi
                  lam[i] = lam[i]  - 2*pi
                end
              end
            end
          end
          for i = 1 : 4
            pts[:,ipts] = [lam[i],theta[i],max(z[i]-Global.Output.RadPrint,0.0)/Global.Output.H/5.0]
#           pts[:,ipts] = [lam[i],theta[i],0.0]
            ipts = ipts + 1
          end
        else
          for i=1:4
            pts[:,ipts] = [x[i,1],x[i,2],x[i,3]]
            ipts = ipts + 1
          end
        end
        ksi0=ksi1
      end
      eta0=eta1
    end
  end
  celltype = VTKCellTypes.VTK_QUAD 

  ConnectivityList=reshape(1:1:4*NF*OrdPrint*OrdPrint,4,NF*OrdPrint*OrdPrint)
  cells = MeshCell[]

  for k in 1 : NF * OrdPrint * OrdPrint
    inds = Vector(1 : 4) .+ 4 * (k -1)
    push!(cells, MeshCell(celltype, inds))
  end

  vtkInter = zeros(Float64,OrdPrint,OrdPrint,1,OrdPoly+1,OrdPoly+1,1)
  dd=2/OrdPrint
  eta0=-1
  for jRef=1:OrdPrint
    ksi0=-1
    eta1=eta0+dd
    for iRef=1:OrdPrint
      ksi1=ksi0+dd
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          vtkInter[iRef,jRef,1,i,j,1] = DG.Lagrange(0.5*(ksi0+ksi1),FE.xwCPU,i)*
                DG.Lagrange(0.5*(eta0+eta1),FE.xwCPU,j)
        end
      end
      ksi0 = ksi1
    end
    eta0 = eta1
  end
  dvtkInter = KernelAbstractions.zeros(backend,FT,size(vtkInter))
  copyto!(dvtkInter,vtkInter)
  RefineMidPoints = zeros(0,0)
  return vtkStruct{FT,
                   typeof(dvtkInter)}(
    dvtkInter,
    cells,
    pts,
    RefineMidPoints,
  )  
end

function vtkSkeleton!(vtkCache,filename, part::Int, nparts::Int, c, FileNumber, cName)
  cells = vtkCache.cells
  pts = vtkCache.pts

  step = FileNumber
  stepS = "$step"
  vtk_filename_noext = pwd()*"/output/VTK/" * filename * stepS
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)
  for iC = 1 : length(cName)
    vtk[cName[iC], VTKCellData()] = c[:,iC]
  end
  outfiles = vtk_save(vtk)
  return nothing
end  

function unstructured_vtkSphere(U,Trans,FE,Metric,Phys,Global, part::Int, nparts::Int;
  Cache=zeros(0,0,0))

  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  NG = FE.NumG
  OrdPoly = FE.OrdPoly 
  OrdPolyZ = FE.OrdPolyZ 
  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  OrdPrint = Global.Output.OrdPrint
  OrdPrintZ = Global.Output.OrdPrintZ
  vtkInter = Global.vtkCache.vtkInter
  cells = Global.vtkCache.cells
  pts = Global.vtkCache.pts
  filename = Global.Output.vtkFileName

  step = Global.Output.vtk
  stepS="$step"
  vtk_filename_noext = filename * stepS
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)

  backend = get_backend(U)				      
  FTB = eltype(U)
  if length(size(U)) == 3
    UR = reshape(U,size(U,1),1,size(U,2),size(U,3))  
  else
    UR = U
  end  
  cCell = KernelAbstractions.zeros(backend,FTB,OrdPrint,OrdPrint,OrdPrintZ,nz,NF)
  cCellCPU = zeros(OrdPrint*OrdPrint*OrdPrintZ*nz*NF) 
  for i=1:length(Global.Output.cNames)
    str = Global.Output.cNames[i]
    if str == "Rho"
      RhoPos = Global.Model.RhoPos
      @views InterpolateGPU!(cCell,UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["rho", VTKCellData()] = cCellCPU
    elseif  str == "u" 
      uPos = Global.Model.uPos
      @views InterpolateGPU!(cCell,UR[:,:,:,uPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["u", VTKCellData()] = cCellCPU
    elseif  str == "Rhou" 
      uPos = Global.Model.uPos
      RhoPos = Global.Model.RhoPos
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,uPos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["u", VTKCellData()] = cCellCPU
    elseif  str == "v" 
      vPos = Global.Model.vPos
      @views InterpolateGPU!(cCell,UR[:,:,:,vPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["v", VTKCellData()] = cCellCPU
    elseif  str == "Thermo" 
      ThPos = Global.Model.ThPos
      @views InterpolateGPU!(cCell,UR[:,:,:,ThPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["Thermo", VTKCellData()] = cCellCPU
    elseif  str == "Rhov" 
      vPos = Global.Model.vPos
      RhoPos = Global.Model.RhoPos
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,vPos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["v", VTKCellData()] = cCellCPU
    elseif  str == "Rhow" 
      wPos = Global.Model.wPos
      RhoPos = Global.Model.RhoPos
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,wPos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["w", VTKCellData()] = cCellCPU
    elseif  str == "wDG" 
      wPos = Global.Model.wPos
      wCell = zeros(OrdPrint*OrdPrint*OrdPrintZ*nz*NF)
      @. @views UR[:,:,FE.BoundaryDoF,wPos] = FTB(0.0)
      @views InterpolateGPU!(cCell,UR[:,:,:,wPos],vtkInter,FE.Glob)
      @views copyto!(wCell,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["wDG", VTKCellData()] = wCell
    elseif  str == "w" 
      uPos = Global.Model.uPos
      vPos = Global.Model.uPos
      wPos = Global.Model.wPos
      wCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateW!(wCell,UR[:,:,:,wPos],UR[:,:,:,uPos],UR[:,:,:,vPos],
        vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz,Metric.dXdxI)
      vtk["w", VTKCellData()] = wCell
    elseif  str == "wB" 
      uPos = Global.Model.uPos
      vPos = Global.Model.uPos
      wPos = Global.Model.wPos
      wCell = zeros(OrdPrint*OrdPrint*nz*NF)
      InterpolateWBGPU!(cCell,UR[:,:,:,uPos],UR[:,:,:,vPos],UR[:,:,:,wPos],vtkInter,Metric.dXdxI,FE.Glob)
#     @views InterpolateWB!(wCell,U[:,:,wPos],U[:,:,uPos],U[:,:,vPos],
#       vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz,Metric.dXdxI)
      @views copyto!(wCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["w", VTKCellData()] = wCell
    elseif  str == "BDG" 
      BPos = Global.Model.ThPos
      @views InterpolateGPU!(cCell,UR[:,:,:,BPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["BDG", VTKCellData()] = cCellCPU
    elseif str == "Th"  
      ThCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views PotT = reshape(Cache.Thermo[:,:,3],size(Cache.Thermo[:,:,3],1),1,
        size(Cache.Thermo[:,:,3],2))
      InterpolateGPU!(cCell,PotT,vtkInter,FE.Glob)
      copyto!(ThCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Th", VTKCellData()] = ThCell
    elseif str == "RhoTh"  
      ThPos = Global.Model.ThPos
      RhoPos = Global.Model.RhoPos
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,ThPos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      @views copyto!(cCellCPU,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["Th", VTKCellData()] = cCellCPU
    elseif str == "ThE"  
      ThECell = zeros(OrdPrint*OrdPrint*nz*NF)  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
      else
        RhoPos = Global.Model.RhoPos
        ThPos = Global.Model.ThPos
        RhoTPos = Global.Model.RhoTPos
        RhoVPos = Global.Model.RhoVPos
        RhoCPos = Global.Model.RhoCPos
#       @views InterpolateThE!(ThECell,U[:,:,ThPos],U[:,:,RhoPos],U[:,:,RhoVPos],U[:,:,RhoCPos],
#         vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz,Global.Phys)
        if RhoVPos > 0
          @views InterpolateThEGPU!(cCell,UR[:,:,:,ThPos],UR[:,:,:,RhoPos],UR[:,:,:,RhoVPos],UR[:,:,:,RhoCPos],
            vtkInter,FE.Glob,Phys)
        elseif RhoTPos >0
        end
        copyto!(ThECell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
        vtk["ThE", VTKCellData()] = ThECell 
      end    
    elseif str == "ThDiff"  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
        @views Pres = Cache.Thermo[:,:,1]  
        RhoPos = Global.Model.RhoPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)  
        @views InterpolateTh!(ThCell,Pres,UR[:,:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz,Global.Phys)
        vtk["Th", VTKCellData()] = ThCell 
      else
        RhoPos = Global.Model.RhoPos
        ThPos = Global.Model.ThPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)
        ThCellBGrd = zeros(OrdPrint*OrdPrint*nz*NF)
        @views Interpolate!(ThCell,UR[:,:,:,ThPos],UR[:,:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz)
        @views Interpolate!(ThCellBGrd,Global.ThetaBGrd,vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz)
        @. ThCell -= ThCellBGrd
        vtk["ThDiff", VTKCellData()] = ThCell
      end
    elseif str == "Pres"   
      @views Pres = Cache.Thermo[:,:,1]  
      pCell = zeros(OrdPrint*OrdPrint*OrdPrintZ*nz*NF)
#     Interpolate!(pCell,Pres,vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz)
      if length(size(Pres)) == 2
        PresR = reshape(Pres,size(Pres,1),1,size(Pres,2))  
      else
        PresR = Pres  
      end    
      InterpolateGPU!(cCell,PresR,vtkInter,FE.Glob)
      @views copyto!(pCell,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["Pres", VTKCellData()] = pCell 
    elseif  str == "Tke" 
      TkePos = Global.Model.TkePos
      RhoPos = Global.Model.RhoPos
      TkeCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,TkePos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      copyto!(TkeCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Tke", VTKCellData()] = TkeCell
    elseif  str == "qV" 
      RhoVPos = Global.Model.RhoVPos
      RhoPos = Global.Model.RhoPos
      qVCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,RhoVPos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      copyto!(qVCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["qV", VTKCellData()] = qVCell
    elseif  str == "qC" 
      RhoCPos = Global.Model.RhoCPos
      RhoPos = Global.Model.RhoPos
      qCCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,RhoCPos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      copyto!(qCCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["qC", VTKCellData()] = qCCell
    elseif  str == "DiffKoeff" 
      DiffKoeff = Cache.KV  
      RhoPos = Global.Model.RhoPos
      DiffKoeffCell = zeros(OrdPrint*OrdPrint*OrdPrintZ*nz*NF)
      if length(size(U)) == 3
        DiffKoeffR = reshape(DiffKoeff,size(DiffKoeff,1),1,size(DiffKoeff,2))  
      else
        DiffKoeffR = DiffKoeff  
      end  
      @views InterpolateRhoGPU!(cCell,DiffKoeffR,UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      copyto!(DiffKoeffCell,reshape(cCell,OrdPrint*OrdPrint*OrdPrintZ*nz*NF))
      vtk["DiffKoeff", VTKCellData()] = DiffKoeffCell
    elseif  str == "Tr1" 
      Tr1Pos = Global.Model.NumV + 1
      RhoPos = Global.Model.RhoPos
      Tr1Cell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,Tr1Pos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      copyto!(Tr1Cell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Tr1", VTKCellData()] = Tr1Cell
    elseif  str == "Tr2" 
      Tr2Pos = Global.Model.NumV + 2
      Tr2Cell = zeros(OrdPrint*OrdPrint*nz*NF)
      RhoPos = Global.Model.RhoPos
#     @views Interpolate!(Tr2Cell,U[:,:,Tr2Pos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz)
      @views InterpolateRhoGPU!(cCell,UR[:,:,:,Tr2Pos],UR[:,:,:,RhoPos],vtkInter,FE.Glob)
      copyto!(Tr2Cell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Tr2", VTKCellData()] = Tr2Cell
    elseif str == "Vort"
      uPos = Global.Model.uPos
      vPos = Global.Model.vPos
      VortCell = zeros(OrdPrint*OrdPrint*nz*NF)
#     @views InterpolateVort!(VortCell,U[:,:,uPos:vPos],vtkInter,OrdPrint,FE,Metric,Cache,Global)
      InterpolateVortGPU!(cCell,U,vtkInter,FE,Metric)
      copyto!(VortCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Vort", VTKCellData()] = VortCell
    end   
  end   
  outfiles=vtk_save(vtk)
  Global.Output.vtk = Global.Output.vtk + 1
  return outfiles::Vector{String}
end

function InterpolateWB!(cCell,w,u,v,Inter,OrdPoly,OrdPrint,Glob,NF,nz,dXdxI)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      if iz == 1
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            w0 = -(u[iz,Glob[i,j,iF]] * dXdxI[i,j,1,1,3,1,iF] +
              v[iz,Glob[i,j,iF]] * dXdxI[i,j,1,1,3,2,iF]) / dXdxI[i,j,1,1,3,3,iF]
            @views @. cc = cc + 0.5 * Inter[:,:,i,j]*(w0 + w[2,Glob[i,j,iF]])
          end
        end
      else
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            @views @. cc = cc + 0.5 * Inter[:,:,i,j]*(w[iz-1,Glob[i,j,iF]] + w[iz,Glob[i,j,iF]])
          end
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end



function InterpolateW!(cCell,w,u,v,Inter,OrdPoly,OrdPrint,Glob,NF,nz,dXdxI)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      if iz == 1
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            w0 = 0.0  
            @views @. cc = cc + 0.5 * Inter[:,:,i,j]*(w0 + w[2,Glob[i,j,iF]])
          end
        end  
      else    
        for j=1:OrdPoly+1
          for i=1:OrdPoly+1
            @views @. cc = cc + 0.5 * Inter[:,:,i,j]*(w[iz-1,Glob[i,j,iF]] + w[iz,Glob[i,j,iF]])
          end
        end
      end  
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function Interpolate!(cCell,c,Inter,OrdPoly,OrdPrint,Glob,NF,nz)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      ID = 0
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          ID += 1  
          @views @. cc = cc + Inter[:,:,i,j]*c[iz,Glob[ID,iF]]
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolateClipp!(cCell,c,Rho,Inter,OrdPoly,OrdPrint,Glob,NF,nz)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      minC = 1.e40
      maxC = -1.e40
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          minC = min(minC,c[iz,Glob[i,j,iF]] / Rho[iz,Glob[i,j,iF]])  
          maxC = max(maxC,c[iz,Glob[i,j,iF]] / Rho[iz,Glob[i,j,iF]])  
          @views @. cc = cc + Inter[:,:,i,j]*c[iz,Glob[i,j,iF]] / Rho[iz,Glob[i,j,iF]]
        end
      end
      @. cc = min(cc,maxC)
      @. cc = max(cc,minC)
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolateFE!(cCell,cFE,Inter,OrdPoly,OrdPrint,Glob,NF,nz)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    @. cc = 0.0
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        @views @. cc = cc + Inter[:,:,i,j]*cFE[i,j,iF]
      end
    end
    @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
    icCell = icCell + OrdPrint*OrdPrint
  end
end


function Interpolate!(cCell,c,Rho,Inter,OrdPoly,OrdPrint,Glob,NF,nz)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          @views @. cc = cc + Inter[:,:,i,j]*c[iz,Glob[i,j,iF]] / Rho[iz,Glob[i,j,iF]]
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolateTh!(cCell,Pres,Rho,Inter,OrdPoly,OrdPrint,Glob,NF,nz,Phys)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          cLoc = Pres[iz,Glob[i,j,iF]] / (Rho[iz,Glob[i,j,iF]] * Phys.Rd) * 
            (Phys.p0 / Pres[iz,Glob[i,j,iF]])^Phys.kappa   
          @views @. cc = cc + Inter[:,:,i,j]*cLoc
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolateThE!(cCell,RhoTh,Rho,RhoV,RhoC,Inter,OrdPoly,OrdPrint,Glob,NF,nz,Phys)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          ind = Glob[i,j,iF]  
          cLoc = fThE(Rho[iz,ind],RhoV[iz,ind],RhoC[iz,ind],RhoTh[iz,ind],Phys)
          @views @. cc = cc + Inter[:,:,i,j]*cLoc
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolateVort!(cCell,U,Inter,OrdPrint,FE,Metric,Cache,Global)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  v1FE = Cache.v1FE
  v2FE = Cache.v2FE
  VortFE = Cache.KE
  OP = FE.OrdPoly + 1
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  @inbounds for iF=1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = FE.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          v1FE[iP,jP,iz] = U[iz,ind,1]
          v2FE[iP,jP,iz] = U[iz,ind,2]
        end
      end
    end
    @views Rot!(VortFE,v1FE,v2FE,FE,Metric.dXdxI[:,:,:,:,:,:,iF],
      Metric.J[:,:,:,:,iF],Global.ThreadCache)

    @inbounds for iz=1:nz
      @. cc = 0.0
      @inbounds for j=1:OP
        @inbounds for i=1:OP
          ind = FE.Glob[i,j,iF]  
          @views @. cc = cc + Inter[:,:,i,j]*VortFE[i,j,iz]
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolatePressure!(cCell,U,Inter,OrdPrint,FE,Metric,Cache,Global)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  v1FE = Cache.v1FE
  v2FE = Cache.v2FE
  wFE = Cache.wFE
  wCFE = Cache.wCFE
  KE = Cache.KE
  OP = FE.OrdPoly + 1
  nz = Global.Grid.nz
  NF = Global.Grid.NumFaces
  zP = Metric.zP
  if Global.Model.Thermo == "Energy"
    @inbounds for iF=1:NF
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = FE.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            v1FE[iP,jP,iz] = U[iz,ind,1]
            v2FE[iP,jP,iz] = U[iz,ind,2]
            wFE[iP,jP,iz+1] = U[iz,ind,3]
          end
        end
      end
      KineticEnergy!(KE,v1FE,v2FE,wFE,FE,Global,iF)  
      for iz=1:nz
        @. cc = 0.0
        for j=1:OP
          for i=1:OP
            ind = FE.Glob[i,j,iF]  
            p = Pressure(U[iz,ind,:],KE[i,j,iz],zP[iz,ind],Metric,Global)  
            @views @. cc = cc + Inter[:,:,i,j]*p
          end
        end
        @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
        icCell = icCell + OrdPrint*OrdPrint
      end
    end
  end
end

function unstructured_vtkPartition(vtkGrid, NF, part::Int, nparts::Int)

  nz = 1
  OrdPrint = 1
  vtkInter = vtkGrid.vtkInter
  cells = vtkGrid.cells
  pts = vtkGrid.pts
  filename = "Partition"

  vtk_filename_noext = filename
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)
  PartitionCell = zeros(NF)
  PartitionCell .= part
  vtk["Part", VTKCellData()] = PartitionCell
  outfiles=vtk_save(vtk)
  return outfiles::Vector{String}
end

function unstructured_vtk2Dim(c,vtkGrid, NF, FE,  part::Int, nparts::Int, FileNumber, FileName, cName)
  nz = 1
  OrdPrint = FE.OrdPoly
  OrdPoly = FE.OrdPoly
  vtkInter = vtkGrid.vtkInter
  cells = vtkGrid.cells
  pts = vtkGrid.pts

  step = FileNumber
  stepS = "$step"
  vtk_filename_noext = pwd()*"/output/VTK/" * FileName * stepS
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)
  backend = get_backend(c)
  FTB = eltype(c)
  cCell = KernelAbstractions.zeros(backend,FTB,OrdPrint,OrdPrint,NF)
  HeightCell = zeros(OrdPrint*OrdPrint*NF,length(cName)) 
  iC = 1
  @views InterpolateFEDim2GPU!(cCell,c[:,:,:,iC],vtkInter,FE.Glob)
  @views copyto!(HeightCell[:,iC],cCell)
  @views vtk[cName[iC], VTKCellData()] = HeightCell[:,iC]
  for iC = 2 : length(cName)
    @views InterpolateFEDim2GPU!(cCell,c[:,:,:,iC]./c[:,:,:,1],vtkInter,FE.Glob)
    @views copyto!(HeightCell[:,iC],cCell)
    @views vtk[cName[iC], VTKCellData()] = HeightCell[:,iC]
  end
  outfiles=vtk_save(vtk)
  return outfiles::Vector{String}
end  

function unstructured_vtkOrography(Height,vtkGrid, NF, FE,  part::Int, nparts::Int)
  nz = 1
  OrdPrint = FE.OrdPoly
  OrdPoly = FE.OrdPoly
  vtkInter = vtkGrid.vtkInter
  cells = vtkGrid.cells
  pts = vtkGrid.pts
  filename = "Orography"

  vtk_filename_noext = filename
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)
  HeightCell = zeros(OrdPrint*OrdPrint*NF) 
  backend = get_backend(Height)
  FTB = eltype(Height)
  cCell = KernelAbstractions.zeros(backend,FTB,OrdPrint,OrdPrint,NF)
  #@views InterpolateFE!(HeightCell,Height,vtkInter,OrdPoly,OrdPrint,FE.Glob,NF,nz)
  InterpolateFEDim2GPU!(cCell,Height,vtkInter,FE.Glob)
  copyto!(HeightCell,cCell)
  vtk["Height", VTKCellData()] = HeightCell
  outfiles=vtk_save(vtk)
  return outfiles::Vector{String}
end  

function Bilinear(ksi1,ksi2,P1,P2,P3,P4,Form,R)
  x = 0.25 * ((1 - ksi1) * (1 - ksi2) * P1.x +
    (1 + ksi1) * (1 - ksi2) * P2.x +
    (1 + ksi1) * (1 + ksi2) * P3.x +
    (1 - ksi1) * (1 + ksi2) * P4.x) 
  y = 0.25 * ((1 - ksi1) * (1 - ksi2) * P1.y +
    (1 + ksi1) * (1 - ksi2) * P2.y +
    (1 + ksi1) * (1 + ksi2) * P3.y +
    (1 - ksi1) * (1 + ksi2) * P4.y) 
  z = 0.25 * ((1 - ksi1) * (1 - ksi2) * P1.z +
    (1 + ksi1) * (1 - ksi2) * P2.z +
    (1 + ksi1) * (1 + ksi2) * P3.z +
    (1 - ksi1) * (1 + ksi2) * P4.z) 
  if Form == "Sphere"
    r = sqrt(x^2 + y^2 + z^2)
    x = x / r * R
    y = y / r * R
    z = z / r * R
  end
  return x, y, z
end  
function Linear(ksi1,ksi2,P1,P2,P3,Form,R)
  x = 0.5 * ((-ksi1 - ksi2) * P1.x +
    (1 + ksi1) * P2.x +
    (1 + ksi2) * P3.x)
  y = 0.5 * ((-ksi1 - ksi2) * P1.y +
    (1 + ksi1) * P2.y +
    (1 + ksi2) * P3.y)
  z = 0.5 * ((-ksi1 - ksi2) * P1.z +
    (1 + ksi1) * P2.z +
    (1 + ksi2) * P3.z)
  if Form == "Sphere"
    r = sqrt(x^2 + y^2 + z^2)
    x = x / r * R
    y = y / r * R
    z = z / r * R
  end
  return x, y, z
end
