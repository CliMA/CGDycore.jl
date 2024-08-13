mutable struct vtkStruct{FT<:AbstractFloat,
                         AT4<:AbstractArray}
  vtkInter::AT4
  cells
  pts::Array{Float64,2}
end
function vtkStruct{FT}(backend) where FT<:AbstractFloat
  vtkInter = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  pts = Array{Float64,2}(undef,0,0)
  cells = MeshCell[]
  return vtkStruct{FT,
                   typeof(vtkInter)}(
  vtkInter,
  cells,
  pts,
  )
end

function vtkStruct{FT}(backend,Grid) where FT<:AbstractFloat

  vtkInter = KernelAbstractions.zeros(backend,FT,0,0,0,0)
  celltype = VTKCellTypes.VTK_POLYGON
  cells = MeshCell[]
  NumNodes = 0

  for iF in 1 : Grid.NumFaces
    inds = Vector(1 : length(Grid.Faces[iF].N)) .+ NumNodes
    push!(cells, MeshCell(celltype, inds))
    NumNodes += length(Grid.Faces[iF].N)
  end

  pts = Array{Float64,2}(undef,3,NumNodes)
  NumNodes = 0
  Flat = true
  dTol = 2*pi / 30
  if Flat
    for iF in 1 : Grid.NumFaces
      lam = zeros(length(Grid.Faces[iF].N))
      theta = zeros(length(Grid.Faces[iF].N))
      NumNodesLoc = 0
      for iN in Grid.Faces[iF].N  
        NumNodesLoc += 1  
        (lam[NumNodesLoc],theta[NumNodesLoc],z) = Grids.cart2sphere(Grid.Nodes[iN].P.x,Grid.Nodes[iN].P.y,Grid.Nodes[iN].P.z)
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
    end
  else    
    for iF in 1 : Grid.NumFaces
      for iN in Grid.Faces[iF].N  
        NumNodes += 1  
        pts[1,NumNodes] = Grid.Nodes[iN].P.x  
        pts[2,NumNodes] = Grid.Nodes[iN].P.y  
        pts[3,NumNodes] = Grid.Nodes[iN].P.z  
      end  
    end
  end  

  return vtkStruct{FT,
                   typeof(vtkInter)}(
  vtkInter,
  cells,
  pts,
  )
end                   


function vtkStruct{FT}(backend,OrdPrint::Int,Trans,CG,Metric,Global) where FT<:AbstractFloat
  OrdPoly = CG.OrdPoly
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
          @views Trans(x[1,:],ksi0,eta0, -1.0,X,CG,Global)
          @views Trans(x[2,:],ksi1,eta0, -1.0,X,CG,Global)
          @views Trans(x[3,:],ksi1,eta1, -1.0,X,CG,Global)
          @views Trans(x[4,:],ksi0,eta1, -1.0,X,CG,Global)
          @views Trans(x[5,:],ksi0,eta0, 1.0,X,CG,Global)
          @views Trans(x[6,:],ksi1,eta0, 1.0,X,CG,Global)
          @views Trans(x[7,:],ksi1,eta1, 1.0,X,CG,Global)
          @views Trans(x[8,:],ksi0,eta1, 1.0,X,CG,Global)
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
          vtkInter[iRef,jRef,i,j] = vtkInter[iRef,jRef,i,j] + DG.Lagrange(0.5*(ksi0+ksi1),CG.xwCPU,i)*
                DG.Lagrange(0.5*(eta0+eta1),CG.xwCPU,j)
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

function vtkInit2D(OrdPrint::Int,Trans,CG,Metric,Global)
  OrdPoly = CG.OrdPoly
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
  X = zeros(FT,OrdPoly+1,OrdPoly+1,2,3)
  for iF = 1 : NF
    @views copyto!(X,reshape(Metric.X[:,:,:,1,iF],OrdPoly+1,OrdPoly+1,2,3))
    dd = 2 / OrdPrint
    eta0 = -1
    for jRef = 1 : OrdPrint
      ksi0 = -1
      eta1 = eta0 + dd
      for iRef = 1 : OrdPrint
        ksi1 = ksi0 + dd
        @views Trans(x[1,:],ksi0,eta0, -1.0,X,CG,Global)
        @views Trans(x[2,:],ksi1,eta0, -1.0,X,CG,Global)
        @views Trans(x[3,:],ksi1,eta1, -1.0,X,CG,Global)
        @views Trans(x[4,:],ksi0,eta1, -1.0,X,CG,Global)
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
          vtkInter[iRef,jRef,i,j] = vtkInter[iRef,jRef,i,j] + DG.Lagrange(0.5*(ksi0+ksi1),CG.xw,i)*
                DG.Lagrange(0.5*(eta0+eta1),CG.xw,j)
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

function vtkSkeleton!(vtkCache,filename, part::Int, nparts::Int, c, FileNumber)
  cells = vtkCache.cells
  pts = vtkCache.pts

  step = FileNumber
  stepS = "$step"
  vtk_filename_noext = filename * stepS;
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)
  cName=["Height","uC","vC","wC","uS","vS"]
  for iC = 1 : size(c,2)
    vtk[cName[iC], VTKCellData()] = c[:,iC]
  end
  outfiles = vtk_save(vtk);
  return nothing
end  

function unstructured_vtkSphere(U,Trans,CG,Metric,Cache,Phys,Global, part::Int, nparts::Int)

  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  NG = CG.NumG
  OrdPoly = CG.OrdPoly 
  NumV = Global.Model.NumV
  NumTr = Global.Model.NumTr
  OrdPrint = Global.Output.OrdPrint
  vtkInter = Global.vtkCache.vtkInter
  cells = Global.vtkCache.cells
  pts = Global.vtkCache.pts
  filename = Global.Output.vtkFileName

  step = Global.Output.vtk
  stepS="$step"
  vtk_filename_noext = filename * stepS;
  vtk = pvtk_grid(vtk_filename_noext, pts, cells; compress=3, part = part, nparts = nparts)

  backend = get_backend(U)				      
  FTB = eltype(U)
  cCell = KernelAbstractions.zeros(backend,FTB,OrdPrint,OrdPrint,nz,NF)
  for i=1:length(Global.Output.cNames)
    str = Global.Output.cNames[i]
    if str == "Rho"
      RhoPos = Global.Model.RhoPos
      RhoCell = zeros(OrdPrint*OrdPrint*nz*NF) 
#     @views Interpolate!(RhoCell,Rho,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      @views InterpolateGPU!(cCell,U[:,:,RhoPos],vtkInter,CG.Glob)
      @views copyto!(RhoCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["rho", VTKCellData()] = RhoCell
    elseif  str == "u" 
      uPos = Global.Model.uPos
      uCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateGPU!(cCell,U[:,:,uPos],vtkInter,CG.Glob)
      @views copyto!(uCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["u", VTKCellData()] = uCell
    elseif  str == "Rhou" 
      uPos = Global.Model.uPos
      @views copyto!(temp,U[:,:,uPos])
      uCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views Interpolate!(uCell,temp,Rho,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["u", VTKCellData()] = uCell
    elseif  str == "v" 
      vPos = Global.Model.vPos
      vCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateGPU!(cCell,U[:,:,vPos],vtkInter,CG.Glob)
      @views copyto!(vCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["v", VTKCellData()] = vCell
    elseif  str == "Rhov" 
      vPos = Global.Model.vPos
      RhoPos = Global.Model.RhoPos
      vCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views Interpolate!(vCell,U[:,:,vPos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["v", VTKCellData()] = vCell
    elseif  str == "w" 
      uPos = Global.Model.uPos
      vPos = Global.Model.uPos
      wPos = Global.Model.wPos
      wCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateW!(wCell,U[:,:,wPos],U[:,:,uPos],U[:,:,vPos],
        vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Metric.dXdxI)
      vtk["w", VTKCellData()] = wCell
    elseif  str == "wB" 
      uPos = Global.Model.uPos
      vPos = Global.Model.uPos
      wPos = Global.Model.wPos
      wCell = zeros(OrdPrint*OrdPrint*nz*NF)
      InterpolateWBGPU!(cCell,U[:,:,uPos],U[:,:,vPos],U[:,:,wPos],vtkInter,Metric.dXdxI,CG.Glob)
#     @views InterpolateWB!(wCell,U[:,:,wPos],U[:,:,uPos],U[:,:,vPos],
#       vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Metric.dXdxI)
      @views copyto!(wCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["w", VTKCellData()] = wCell
    elseif str == "Th"  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
        @views Pres = Cache.AuxG[:,:,1]  
        RhoPos = Global.Model.RhoPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)  
        @views InterpolateTh!(ThCell,Pres,U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Phys)
        vtk["Th", VTKCellData()] = ThCell 
      else
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)
        @views PotT = Cache.AuxG[:,:,3]  
        InterpolateGPU!(cCell,PotT,vtkInter,CG.Glob)
        copyto!(ThCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
        vtk["Th", VTKCellData()] = ThCell
      end
    elseif str == "ThE"  
      ThECell = zeros(OrdPrint*OrdPrint*nz*NF)  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
      else
        RhoPos = Global.Model.RhoPos
        ThPos = Global.Model.ThPos
        RhoVPos = Global.Model.RhoVPos
        RhoCPos = Global.Model.RhoCPos
#       @views InterpolateThE!(ThECell,U[:,:,ThPos],U[:,:,RhoPos],U[:,:,RhoVPos],U[:,:,RhoCPos],
#         vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Phys)
        @views InterpolateThEGPU!(cCell,U[:,:,ThPos],U[:,:,RhoPos],U[:,:,RhoVPos],U[:,:,RhoCPos],
          vtkInter,CG.Glob,Phys)
        copyto!(ThECell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
        vtk["ThE", VTKCellData()] = ThECell 
      end    
    elseif str == "ThDiff"  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
        @views Pres = Cache.AuxG[:,:,1]  
        RhoPos = Global.Model.RhoPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)  
        @views InterpolateTh!(ThCell,Pres,U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Phys)
        vtk["Th", VTKCellData()] = ThCell 
      else
        RhoPos = Global.Model.RhoPos
        ThPos = Global.Model.ThPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)
        ThCellBGrd = zeros(OrdPrint*OrdPrint*nz*NF)
        @views Interpolate!(ThCell,U[:,:,ThPos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
        @views Interpolate!(ThCellBGrd,Global.ThetaBGrd,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
        @. ThCell -= ThCellBGrd
        vtk["ThDiff", VTKCellData()] = ThCell
      end
    elseif str == "Pres"   
      @views Pres = Cache.AuxG[:,:,1]  
      pCell = zeros(OrdPrint*OrdPrint*nz*NF)
#     Interpolate!(pCell,Pres,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      InterpolateGPU!(cCell,Pres,vtkInter,CG.Glob)
      @views copyto!(pCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["p", VTKCellData()] = pCell 
    elseif  str == "Tke" 
      TkePos = Global.Model.TkePos
      RhoPos = Global.Model.RhoPos
      TkeCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,U[:,:,TkePos],U[:,:,RhoPos],vtkInter,CG.Glob)
      copyto!(TkeCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Tke", VTKCellData()] = TkeCell
    elseif  str == "DiffKoeff" 
      DiffKoeff = Cache.KV  
      RhoPos = Global.Model.RhoPos
      DiffKoeffCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,DiffKoeff,U[:,:,RhoPos],vtkInter,CG.Glob)
      copyto!(DiffKoeffCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["DiffKoeff", VTKCellData()] = DiffKoeffCell
    elseif  str == "Tr1" 
      Tr1Pos = Global.Model.NumV + 1
      RhoPos = Global.Model.RhoPos
      Tr1Cell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateRhoGPU!(cCell,U[:,:,Tr1Pos],U[:,:,RhoPos],vtkInter,CG.Glob)
      copyto!(Tr1Cell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Tr1", VTKCellData()] = Tr1Cell
    elseif  str == "Tr2" 
      Tr2Pos = Global.Model.NumV + 2
      Tr2Cell = zeros(OrdPrint*OrdPrint*nz*NF)
      RhoPos = Global.Model.RhoPos
#     @views Interpolate!(Tr2Cell,U[:,:,Tr2Pos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      @views InterpolateRhoGPU!(cCell,U[:,:,Tr2Pos],U[:,:,RhoPos],vtkInter,CG.Glob)
      copyto!(Tr2Cell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Tr2", VTKCellData()] = Tr2Cell
    elseif str == "Vort"
      uPos = Global.Model.uPos
      vPos = Global.Model.vPos
      VortCell = zeros(OrdPrint*OrdPrint*nz*NF)
#     @views InterpolateVort!(VortCell,U[:,:,uPos:vPos],vtkInter,OrdPrint,CG,Metric,Cache,Global)
      InterpolateVortGPU!(cCell,U,vtkInter,CG,Metric)
      copyto!(VortCell,reshape(cCell,OrdPrint*OrdPrint*nz*NF))
      vtk["Vort", VTKCellData()] = VortCell
    end   
  end   
  outfiles=vtk_save(vtk);
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

function InterpolateCG!(cCell,cCG,Inter,OrdPoly,OrdPrint,Glob,NF,nz)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  for iF=1:NF
    @. cc = 0.0
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        @views @. cc = cc + Inter[:,:,i,j]*cCG[i,j,iF]
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

function InterpolateVort!(cCell,U,Inter,OrdPrint,CG,Metric,Cache,Global)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  v1CG = Cache.v1CG
  v2CG = Cache.v2CG
  VortCG = Cache.KE
  OP = CG.OrdPoly + 1
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  @inbounds for iF=1:NF
    @inbounds for jP=1:OP
      @inbounds for iP=1:OP
        ind = CG.Glob[iP,jP,iF]
        @inbounds for iz=1:nz
          v1CG[iP,jP,iz] = U[iz,ind,1]
          v2CG[iP,jP,iz] = U[iz,ind,2]
        end
      end
    end
    @views Rot!(VortCG,v1CG,v2CG,CG,Metric.dXdxI[:,:,:,:,:,:,iF],
      Metric.J[:,:,:,:,iF],Global.ThreadCache)

    @inbounds for iz=1:nz
      @. cc = 0.0
      @inbounds for j=1:OP
        @inbounds for i=1:OP
          ind = CG.Glob[i,j,iF]  
          @views @. cc = cc + Inter[:,:,i,j]*VortCG[i,j,iz]
        end
      end
      @views cCell[icCell:icCell+OrdPrint*OrdPrint-1] = reshape(cc,OrdPrint*OrdPrint)
      icCell = icCell + OrdPrint*OrdPrint
    end
  end
end

function InterpolatePressure!(cCell,U,Inter,OrdPrint,CG,Metric,Cache,Global)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  v1CG = Cache.v1CG
  v2CG = Cache.v2CG
  wCG = Cache.wCG
  wCCG = Cache.wCCG
  KE = Cache.KE
  OP = CG.OrdPoly + 1
  nz = Global.Grid.nz
  NF = Global.Grid.NumFaces
  zP = Metric.zP
  if Global.Model.Thermo == "Energy"
    @inbounds for iF=1:NF
      @inbounds for jP=1:OP
        @inbounds for iP=1:OP
          ind = CG.Glob[iP,jP,iF]
          @inbounds for iz=1:nz
            v1CG[iP,jP,iz] = U[iz,ind,1]
            v2CG[iP,jP,iz] = U[iz,ind,2]
            wCG[iP,jP,iz+1] = U[iz,ind,3]
          end
        end
      end
      KineticEnergy!(KE,v1CG,v2CG,wCG,CG,Global,iF)  
      for iz=1:nz
        @. cc = 0.0
        for j=1:OP
          for i=1:OP
            ind = CG.Glob[i,j,iF]  
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
  outfiles=vtk_save(vtk);
  return outfiles::Vector{String}
end

function unstructured_vtkOrography(Height,vtkGrid, NF, CG,  part::Int, nparts::Int)
  nz = 1
  OrdPrint = CG.OrdPoly
  OrdPoly = CG.OrdPoly
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
  #@views InterpolateCG!(HeightCell,Height,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
  InterpolateCGDim2GPU!(cCell,Height,vtkInter,CG.Glob)
  copyto!(HeightCell,cCell)
  vtk["Height", VTKCellData()] = HeightCell
  outfiles=vtk_save(vtk);
  return outfiles::Vector{String}
end  
