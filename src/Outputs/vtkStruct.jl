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
      if Grid.Form == Grids.SphericalGrid()
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
      RefinePoints, RefineMidPoints = Grids.PrintPoints(Refine,Grids.Quad())
      NumRefine = size(RefinePoints,1)
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
          if Grid.Form == Grids.SphericalGrid() 
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
              for i in 1 : 4
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
      RefinePoints, RefineMidPoints = Grids.PrintPoints(Refine,Grids.Tri())
      NumRefine = size(RefinePoints,1)
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
          if Grid.Form == Grids.SphericalGrid() 
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
