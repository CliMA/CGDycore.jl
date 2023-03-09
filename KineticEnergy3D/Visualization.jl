mutable struct vtkStruct
  cells
  pts::Array{Float64,2}
  Step::Int64
end
function Plot2DC(c,Fe,xP,zP,Name)
  Nx = size(c,1)
  Nz = size(c,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  cPlot = zeros(Nx*OrdPolyX,Nz*OrdPolyZ)

  for ix = 1 : Nx
    for iz = 1 : Nz  
      cPlot[(ix-1)*OrdPolyX + 1 : ix*OrdPolyX, (iz-1)*OrdPolyZ + 1 : iz*OrdPolyZ] =
        Fe.IntXF2cE*c[ix,iz,:,:]*Fe.IntZC2cE' 
    end
  end  
  fig = Figure()
  ax = Axis(fig[1, 1])
  hm = heatmap!(ax, cPlot)
  Colorbar(fig[1, 2], hm)
  fig
  save(Name*".png", fig)
end  

function Plot2DF(c,Fe,xP,zP,Name)
  Nx = size(c,1)
  Nz = size(c,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  cPlot = zeros(Nx*OrdPolyX,Nz*OrdPolyZ)
  xPlot = zeros(Nx*OrdPolyX+1,Nz*OrdPolyZ+1)
  zPlot = zeros(Nx*OrdPolyX+1,Nz*OrdPolyZ+1)

  for ix = 1 : Nx
    for iz = 1 : Nz
      @views cPlot[(ix-1)*OrdPolyX + 1 : ix*OrdPolyX, (iz-1)*OrdPolyZ + 1 : iz*OrdPolyZ] =
        Fe.IntXF2cE*c[ix,iz,:,:]*Fe.IntZF2cE'
      @views @. xPlot[(ix-1)*OrdPolyX + 1 : (ix*OrdPolyX + 1), (iz-1)*OrdPolyZ + 1 : (iz*OrdPolyZ + 1)] =  
        xP[ix,iz,:,:]
      @views @. zPlot[(ix-1)*OrdPolyX + 1 : (ix*OrdPolyX + 1), (iz-1)*OrdPolyZ + 1 : (iz*OrdPolyZ + 1)] =  
        zP[ix,iz,:,:]
    end
  end
  points = vec([Point2f(xv, zv) for (xv, zv) in zip(xPlot, zPlot)])
  faces = decompose(QuadFace{GLIndex}, Tesselation(Rect(0, 0, 1, 1), size(xPlot)))
  gb_mesh = GeometryBasics.Mesh(meta(points), faces)
  fig, ax, pl = mesh(gb_mesh,  color = cPlot, colormap=:blues)
  wireframe!(ax, gb_mesh)
  save(Name*".png", fig)
end

function vtkPlot3DC(cC,Fe,vtkS,Name)
  Nx = size(cC,1)
  Ny = size(cC,2)
  Nz = size(cC,3)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  Step = vtkS.Step
  stepS="$Step"
  vtk_filename_noext = Name *stepS
  vtk = pvtk_grid(vtk_filename_noext, vtkS.pts, vtkS.cells; compress=3, part = 1, nparts = 1)
  cPlot = zeros(OrdPolyX,OrdPolyY)
  pCell = zeros(Nx*Ny*Nz*OrdPolyX*OrdPolyY)
  iCell = 1
  iCellIncr = OrdPolyX * OrdPolyY
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        @views ComputeOutputC!(cPlot,cC[ix,iy,iz,:,:],Fe)  
        for i = 1 : OrdPolyX   
          for j = 1 : OrdPolyY   
            pCell[iCell]= cPlot[i,j]
            iCell += 1
          end  
        end  
      end
    end
  end  
  vtk = pvtk_grid(vtk_filename_noext, vtkS.pts, vtkS.cells; compress=3, part = 1, nparts = 1)
  vtk[Name, VTKCellData()] = pCell
  outfiles=vtk_save(vtk);
  return outfiles::Vector{String}
end

function vtkPlot3DF(cF,Fe,vtkS,Name)
  Nx = size(cF,1)
  Ny = size(cF,2)
  Nz = size(cF,3) - 1
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY

  Step = vtkS.Step
  stepS="$Step"
  vtk_filename_noext = Name *stepS
  vtk = pvtk_grid(vtk_filename_noext, vtkS.pts, vtkS.cells; compress=3, part = 1, nparts = 1)
  cPlot = zeros(OrdPolyX,OrdPolyY)
  pCell = zeros(Nx*Ny*Nz*OrdPolyX*OrdPolyY)
  iCell = 1
  iCellIncr = OrdPolyX * OrdPolyY
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        @views ComputeOutputC!(cPlot,0.5*(cF[ix,iy,iz,:,:]+cF[ix,iy,iz+1,:,:]),Fe)  
        for i = 1 : OrdPolyX   
          for j = 1 : OrdPolyY   
            pCell[iCell]= cPlot[i,j]
            iCell += 1
          end  
        end  
      end
    end
  end  
  vtk = pvtk_grid(vtk_filename_noext, vtkS.pts, vtkS.cells; compress=3, part = 1, nparts = 1)
  vtk[Name, VTKCellData()] = pCell
  outfiles=vtk_save(vtk);
  return outfiles::Vector{String}
end

function vtkStruct(Fe,xP,yP,zP)
  Nx = size(xP,1)
  Ny = size(xP,2)
  Nz = size(xP,3)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ

  nCells = Nx * Ny * Nz * OrdPolyX * OrdPolyY * OrdPolyZ
  Npts = 8 * nCells
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  for ix = 1 : Nx
    for iy = 1 : Ny
      for iz = 1 : Nz
        for i = 1 : OrdPolyX  
          for j = 1 : OrdPolyY  
            for k = 1 : OrdPolyZ  
              pts[:,ipts] = [xP[ix,iy,iz,i,j,k],yP[ix,iy,iz,i,j,k],zP[ix,iy,iz,i,j,k]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i+1,j,k],yP[ix,iy,iz,i+1,j,k],zP[ix,iy,iz,i+1,j,k]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i+1,j+1,k],yP[ix,iy,iz,i+1,j+1,k],zP[ix,iy,iz,i+1,j+1,k]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i,j+1,k],yP[ix,iy,iz,i,j+1,k],zP[ix,iy,iz,i,j+1,k]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i,j,k],yP[ix,iy,iz,i,j,k],zP[ix,iy,iz,i,j,k+1]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i+1,j,k],yP[ix,iy,iz,i+1,j,k],zP[ix,iy,iz,i+1,j,k+1]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i+1,j+1,k],yP[ix,iy,iz,i+1,j+1,k],zP[ix,iy,iz,i+1,j+1,k+1]]
              ipts = ipts + 1
              pts[:,ipts] = [xP[ix,iy,iz,i,j+1,k],yP[ix,iy,iz,i,j+1,k],zP[ix,iy,iz,i,j+1,k+1]]
              ipts = ipts + 1
            end  
          end  
        end  
      end  
    end  
  end  
  celltype = VTKCellTypes.VTK_HEXAHEDRON
  ConnectivityList=reshape(1:1:8*nCells,8,nCells)
  cells = MeshCell[]

  for k in 1 : nCells
    inds = Vector(1 : 8) .+ 8 * (k -1)
    push!(cells, MeshCell(celltype, inds))
  end
  Step = 1
  return vtkStruct(
    cells,
    pts,
    Step,
  )
end  
