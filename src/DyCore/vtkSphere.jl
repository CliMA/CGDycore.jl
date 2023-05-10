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

function vtkStruct(OrdPrint::Int,Trans,CG,Global)
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

function vtkInit2D(OrdPrint::Int,Trans,CG,Global)
  OrdPoly = CG.OrdPoly
  NF = Global.Grid.NumFaces
  Npts = 4 * NF * OrdPrint * OrdPrint
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  X = zeros(4,3)
  lam=zeros(4,1)
  theta=zeros(4,1)
  z=zeros(4,1)
  if Global.Grid.Form == "Sphere" && Global.Output.Flat
    dTol=2*pi/max(Global.Output.nPanel,1)/4
  end

  for iF = 1 : NF
    dd = 2 / OrdPrint
    eta0 = -1
    for jRef = 1 : OrdPrint
      ksi0 = -1
      eta1 = eta0 + dd
      for iRef = 1 : OrdPrint
        ksi1 = ksi0 + dd
        X[1,:] = Trans(ksi0,eta0, -1.0,Global.Metric.X[:,:,:,:,1,iF],CG,Global)
        X[2,:] = Trans(ksi1,eta0, -1.0,Global.Metric.X[:,:,:,:,1,iF],CG,Global)
        X[3,:] = Trans(ksi1,eta1, -1.0,Global.Metric.X[:,:,:,:,1,iF],CG,Global)
        X[4,:] = Trans(ksi0,eta1, -1.0,Global.Metric.X[:,:,:,:,1,iF],CG,Global)
        if Global.Grid.Form == "Sphere" && Global.Output.Flat
          for i=1:4
            (lam[i],theta[i],z[i]) = cart2sphere(X[i,1],X[i,2],X[i,3])
          end 
          lammin = minimum(lam)
          lammax = maximum(lam)
          if lammin < 0.0 || lammax > 2*pi
#           @show lammin,lammax  
#           stop
          end  
          if abs(lammin - lammax) > 2*pi-dTol
#           @show "vor",lam
            for i = 1 : 4
              if lam[i] < pi
                lam[i] = lam[i] + 2*pi
                if lam[i] > 3*pi
                  lam[i] = lam[i]  - 2*pi
                end
              end
            end
#           @show "nac",lam
          end
          for i = 1 : 4
            pts[:,ipts] = [lam[i],theta[i],max(z[i]-Global.Output.RadPrint,0.0)/Global.Output.H/5.0]
            ipts = ipts + 1
          end
        else
          for i=1:4
            pts[:,ipts] = [X[i,1],X[i,2],X[i,3]]
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

function unstructured_vtkSphere(U,Trans,CG,Global, part::Int, nparts::Int)

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

  for i=1:length(Global.Output.cNames)
    str = Global.Output.cNames[i]
    if str == "Rho"
      RhoPos = Global.Model.RhoPos
      RhoCell = zeros(OrdPrint*OrdPrint*nz*NF) 
      @views Interpolate!(RhoCell,U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["rho", VTKCellData()] = RhoCell
    elseif  str == "u" 
      uPos = Global.Model.uPos
      uCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views Interpolate!(uCell,U[:,:,uPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["u", VTKCellData()] = uCell
    elseif  str == "Rhou" 
      uPos = Global.Model.uPos
      RhoPos = Global.Model.RhoPos
      uCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views Interpolate!(uCell,U[:,:,uPos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["u", VTKCellData()] = uCell
    elseif  str == "v" 
      vPos = Global.Model.vPos
      vCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views Interpolate!(vCell,U[:,:,vPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
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
        vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Metric.dXdxI)
      vtk["w", VTKCellData()] = wCell
    elseif  str == "wB" 
      uPos = Global.Model.uPos
      vPos = Global.Model.uPos
      wPos = Global.Model.wPos
      wCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateWB!(wCell,U[:,:,wPos],U[:,:,uPos],U[:,:,vPos],
        vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Metric.dXdxI)
      vtk["w", VTKCellData()] = wCell
    elseif str == "Th"  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
        RhoPos = Global.Model.RhoPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)  
        @views InterpolateTh!(ThCell,Global.Cache.PresG,U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Phys)
        vtk["Th", VTKCellData()] = ThCell 
      else
        RhoPos = Global.Model.RhoPos
        ThPos = Global.Model.ThPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)
        @views Interpolate!(ThCell,U[:,:,ThPos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
        vtk["Th", VTKCellData()] = ThCell
      end
    elseif str == "ThDiff"  
      if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
        RhoPos = Global.Model.RhoPos
        ThCell = zeros(OrdPrint*OrdPrint*nz*NF)  
        @views InterpolateTh!(ThCell,Global.Cache.PresG,U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz,Global.Phys)
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
      pCell = zeros(OrdPrint*OrdPrint*nz*NF)
      Interpolate!(pCell,Global.Cache.PresG,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["p", VTKCellData()] = pCell 
    elseif  str == "Tr1" 
      Tr1Pos = Global.Model.NumV + 1
      Tr1Cell = zeros(OrdPrint*OrdPrint*nz*NF)
      RhoPos = Global.Model.RhoPos
      @views Interpolate!(Tr1Cell,U[:,:,Tr1Pos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["Tr1", VTKCellData()] = Tr1Cell
    elseif  str == "Tr2" 
      Tr2Pos = Global.Model.NumV + 2
      Tr2Cell = zeros(OrdPrint*OrdPrint*nz*NF)
      RhoPos = Global.Model.RhoPos
      @views Interpolate!(Tr2Cell,U[:,:,Tr2Pos],U[:,:,RhoPos],vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
      vtk["Tr2", VTKCellData()] = Tr2Cell
    elseif str == "Vort"
      uPos = Global.Model.uPos
      vPos = Global.Model.vPos
      VortCell = zeros(OrdPrint*OrdPrint*nz*NF)
      @views InterpolateVort!(VortCell,U[:,:,uPos:vPos],vtkInter,OrdPrint,CG,Global)
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

function InterpolateVort!(cCell,U,Inter,OrdPrint,CG,Global)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  VortCG = Global.Cache.KE
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
    @views Rot!(VortCG,v1CG,v2CG,CG,Global.Metric.dXdxI[:,:,:,:,:,:,iF],
      Global.Metric.J[:,:,:,:,iF],Global.ThreadCache)

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

function InterpolatePressure!(cCell,U,Inter,OrdPrint,CG,Global)
  icCell  = 1
  cc=zeros(OrdPrint,OrdPrint)
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  KE = Global.Cache.KE
  OP = CG.OrdPoly + 1
  nz = Global.Grid.nz
  NF = Global.Grid.NumFaces
  zP = Global.Metric.zP
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
            p = Pressure(U[iz,ind,:],KE[i,j,iz],zP[iz,ind],Global)  
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
  @views InterpolateCG!(HeightCell,Height,vtkInter,OrdPoly,OrdPrint,CG.Glob,NF,nz)
  vtk["Height", VTKCellData()] = HeightCell
  outfiles=vtk_save(vtk);
  return outfiles::Vector{String}
end  
