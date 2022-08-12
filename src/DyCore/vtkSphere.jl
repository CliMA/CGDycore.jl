mutable struct vtkStruct
  vtkInter::Array{Float64,4}
  cells
  pts::Array{Float64,2}
end
function vtkStruct()
  vtkInter = zeros(Float64,0,0,0,0)
  pts = Array{Float64,2}(undef,0,0)
  cells = MeshCell[]
    return vtkStruct(
    vtkInter,
    cells,
    pts,
  )
end

function vtkInit(OrdPrint::Int,Trans,CG,Global)
  OrdPoly = CG.OrdPoly
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  Npts = 8 * NF * nz * OrdPrint * OrdPrint
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  X = zeros(8,3)
  lam=zeros(8,1)
  theta=zeros(8,1)
  z=zeros(8,1)
  if Global.Grid.Form == "Sphere" && Global.Output.Flat
    dTol=2*pi/max(Global.Output.nPanel-1,1)
  end

  for iF = 1 : NF
    for iz = 1 : nz
      dd = 2 / OrdPrint
      eta0 = -1
      for jRef = 1 : OrdPrint
        ksi0 = -1
        eta1 = eta0 + dd
        for iRef = 1 : OrdPrint
          ksi1 = ksi0 + dd
          X[1,:] = Trans(ksi0,eta0, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[2,:] = Trans(ksi1,eta0, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[3,:] = Trans(ksi1,eta1, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[4,:] = Trans(ksi0,eta1, -1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[5,:] = Trans(ksi0,eta0, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[6,:] = Trans(ksi1,eta0, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[7,:] = Trans(ksi1,eta1, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          X[8,:] = Trans(ksi0,eta1, 1.0,Global.Metric.X[:,:,:,:,iz,iF],CG,Global)
          if Global.Grid.Form == "Sphere" && Global.Output.Flat
            for i=1:8
              (lam[i],theta[i],z[i]) = cart2sphere(X[i,1],X[i,2],X[i,3])
            end 
            lammin = minimum(lam)
            lammax = maximum(lam)
            if abs(lammin - lammax) > 2*pi-dTol
              for i = 1 : 8
                if lam[i] < pi
                  lam[i] = lam[i] + 2*pi
                  if lam[i] > 3*pi
                    lam[i] = lam[i]  - 2*pi
                  end
                end
              end
            end
            for i = 1 : 8
              pts[:,ipts] = [lam[i],theta[i],max(z[i]-Global.Output.RadPrint,0.0)/Global.Output.H*3]
              ipts = ipts + 1
            end
          else
            for i=1:8
              pts[:,ipts] = [X[i,1],X[i,2],X[i,3]]
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

  vtkInter = zeros(Float64,OrdPrint,OrdPrint,OrdPoly+1,OrdPoly+1)
  dd=2/OrdPrint;
  eta0=-1;
  for jRef=1:OrdPrint
    ksi0=-1;
    eta1=eta0+dd;
    for iRef=1:OrdPrint
      ksi1=ksi0+dd;
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          vtkInter[iRef,jRef,i,j] = vtkInter[iRef,jRef,i,j] + Lagrange(0.5*(ksi0+ksi1),CG.xw,i)*
                Lagrange(0.5*(eta0+eta1),CG.xw,j)
        end
      end
      ksi0 = ksi1
    end
    eta0 = eta1
  end
  return vtkStruct(
    vtkInter,
    cells,
    pts,
  )  
end

function unstructured_vtkSphere(U,Trans,CG,Global, filename::String, part::Int, nparts::Int)

  @show "unstructured_vtkSphere"
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  NG = CG.NumG
  OrdPoly = CG.OrdPoly 
  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  RhoPos = Global.Model.RhoPos
  OrdPrint = Global.Output.OrdPrint
  vtkInter = Global.vtkCache.vtkInter
  cells = Global.vtkCache.cells
  pts = Global.vtkCache.pts


  @show Global.Output.vtk
  step = Global.Output.vtk
  stepS="$step"
  vtk_filename_noext = pwd()*"/"*filename * stepS;
  @show vtk_filename_noext
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)


  RhoCell = zeros(OrdPrint*OrdPrint*nz*NF) 
  @views Interpolate!(RhoCell,U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
  vtk["Rho", VTKCellData()] = RhoCell

  uPos = Global.Model.uPos
  uCell = zeros(OrdPrint*OrdPrint*nz*NF)
  @views Interpolate!(uCell,U[:,:,uPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
  vtk["u", VTKCellData()] = uCell

  outfiles=vtk_save(vtk);
  Global.Output.vtk = Global.Output.vtk + 1
  return outfiles::Vector{String}
end

function Interpolate!(cCell,c,Inter,OrdPoly,OrdPrint,Glob,NF,nz)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    for iz=1:nz
      @. cc = 0.0
      for j=1:OrdPoly+1
        for i=1:OrdPoly+1
          @views @. cc = cc + Inter[:,:,i,j]*c[iz,Glob[i,j,iF]]
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

