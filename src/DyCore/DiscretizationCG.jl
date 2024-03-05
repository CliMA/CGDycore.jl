function DiscretizationCG(backend,FT,Jacobi,CG::CGQuad,Exchange,Global,zs) 
# Discretization
  Grid = Global.Grid
  OP = CG.OrdPoly+1
  DoF = CG.DoF
  OPZ = CG.OrdPolyZ+1
  nz = Grid.nz
  NF = Grid.NumFaces

  nQuad = OP * OP
  Metric = MetricStruct{FT}(backend,DoF,OPZ,Global.Grid.NumFaces,nz)
  F = zeros(4,3,NF)
  FGPU = KernelAbstractions.zeros(backend,FT,4,3,NF)
  for iF = 1 : NF
    F[1,1,iF] = Grid.Faces[iF].P[1].x  
    F[1,2,iF] = Grid.Faces[iF].P[1].y  
    F[1,3,iF] = Grid.Faces[iF].P[1].z  
    F[2,1,iF] = Grid.Faces[iF].P[2].x  
    F[2,2,iF] = Grid.Faces[iF].P[2].y  
    F[2,3,iF] = Grid.Faces[iF].P[2].z  
    F[3,1,iF] = Grid.Faces[iF].P[3].x  
    F[3,2,iF] = Grid.Faces[iF].P[3].y  
    F[3,3,iF] = Grid.Faces[iF].P[3].z  
    F[4,1,iF] = Grid.Faces[iF].P[4].x  
    F[4,2,iF] = Grid.Faces[iF].P[4].y  
    F[4,3,iF] = Grid.Faces[iF].P[4].z  
  end  
  copyto!(FGPU,F)
  if Global.Grid.Form == "Sphere"
    Grids.JacobiSphere3GPU!(Global.Grid.AdaptGrid,Metric.X,Metric.dXdx,Metric.dXdxI,Metric.J,CG,FGPU,
      Grid.z,zs,Grid.Rad,Global.Model.Equation)
  else
    Grids.JacobiDG3GPU!(Metric.X,Metric.dXdxI,Metric.J,CG,FGPU,Grid.z,zs)
  end  

  MassCGGPU!(CG,Metric.J,CG.Glob,Exchange,Global)

  Metric.dz = KernelAbstractions.zeros(backend,FT,nz,CG.NumG)
  Metric.zP = KernelAbstractions.zeros(backend,FT,nz,CG.NumG)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU  
  NzG = min(div(NumberThreadGPU,DoF),nz)  
  group = (DoF,NzG,1)  
  ndrange = (DoF,nz,NF)
  if Global.Grid.Form == "Sphere"
    Metric.lat = KernelAbstractions.zeros(backend,FT,CG.NumG)
    KGridSizeSphereKernel! = GridSizeSphereKernel!(backend,group)
    Rad = Global.Grid.Rad
    KGridSizeSphereKernel!(Metric.lat,Metric.zP,Metric.dz,Metric.X,CG.Glob,
      Rad,ndrange=ndrange)
  else
    Metric.lat = KernelAbstractions.zeros(backend,FT,0)
    KGridSizeCartKernel! = GridSizeCartKernel!(backend,group)
    KGridSizeCartKernel!(Metric.zP,Metric.dz,Metric.X,CG.Glob,ndrange=ndrange)
  end    
  NFG = min(div(NumberThreadGPU,DoF),NF)  
  group = (DoF, NFG)
  ndrange = (DoF, NF)
  KSurfaceNormalKernel! = SurfaceNormalKernel!(backend,group)
  KSurfaceNormalKernel!(Metric.FS,Metric.nS,Metric.dXdxI,ndrange=ndrange)

  if Global.Model.HorLimit
    Metric.JC = KernelAbstractions.zeros(backend,FT,size(Metric.J,1),size(Metric.J,3),size(Metric.J,4))  
    Metric.JCW = KernelAbstractions.zeros(backend,FT,size(Metric.J,1),size(Metric.J,3),size(Metric.J,4))  
#   NFG = min(div(NumberThreadGPU,),NF)  
    group = (nz, 1)
    ndrange = (nz, NF)
    KCenterJacobiansKernel! = CenterJacobiansKernel!(backend,group)
    KCenterJacobiansKernel!(CG.OrdPoly+1,Metric.JC,Metric.JCW,Metric.J,CG.w,ndrange=ndrange)
  end    

  return (CG,Metric)
end

function DiscretizationCG(backend,FT,Jacobi,CG,Exchange,Global)
  DiscretizationCG(backend,FT,Jacobi,CG,Exchange,Global,
  KernelAbstractions.zeros(backend,FT,CG.OrdPoly+1,CG.OrdPoly+1,Global.Grid.NumFaces))
end  

@kernel inbounds = true function SurfaceNormalKernel!(FS,nS,@Const(dXdxI))

  ID,IF = @index(Global, NTuple)

  NF = @uniform @ndrange()[2]

  if IF <= NF
    FS[ID,IF] = sqrt(dXdxI[3,1,1,ID,1,IF] * dXdxI[3,1,1,ID,1,IF] +
      dXdxI[3,2,1,ID,1,IF] * dXdxI[3,2,1,ID,1,IF] + 
      dXdxI[3,3,1,ID,1,IF] * dXdxI[3,3,1,ID,1,IF])
    nS[ID,1,IF] = dXdxI[3,1,1,ID,1,IF] / FS[ID,IF]
    nS[ID,2,IF] = dXdxI[3,2,1,ID,1,IF] / FS[ID,IF]
    nS[ID,3,IF] = dXdxI[3,3,1,ID,1,IF] / FS[ID,IF]
  end
end

@kernel inbounds = true function CenterJacobiansKernel!(N,JC,JCW,@Const(J),@Const(w))
  Iz,IF = @index(Global, NTuple)

  FT = eltype(JC)
  Nz = @uniform @ndrange()[1]
  if Iz <= Nz
    @views sumJ = sum(J[:,:,Iz,IF])
    for j = 1 : N
      for i = 1 : N
        ID = i + (j - 1) * N  
        JC[ID,Iz,IF] = J[ID,1,Iz,IF] + J[ID,2,Iz,IF] 
        JCW[ID,Iz,IF] = JC[ID,Iz,IF] * w[i] * w[j] / sumJ
      end
    end  
  end  
end

@kernel inbounds = true function GridSizeSphereKernel!(lat,zP,dz,@Const(X),@Const(Glob),Rad)

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]


  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]
    x = eltype(X)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(X)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(X)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    r = sqrt(x * x + y * y + z * z)
    if Iz == 1
      lat[ind] = asin(z / r)
    end  
    zP[Iz,ind] = max(r-Rad, eltype(X)(0))
    x = X[ID,1,1,Iz,IF]
    y = X[ID,1,2,Iz,IF]
    z = X[ID,1,3,Iz,IF]
    r1 = sqrt(x * x + y * y + z * z)
    x = X[ID,2,1,Iz,IF]
    y = X[ID,2,2,Iz,IF]
    z = X[ID,2,3,Iz,IF]
    r2 = sqrt(x * x + y * y + z * z)
    dz[Iz,ind] =  r2 - r1
  end
end

@kernel inbounds = true function GridSizeCartKernel!(zP,dz,@Const(X),@Const(Glob))

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]
  
  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]  
    zP[Iz,ind] = eltype(X)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF]) 
    dz[Iz,ind] = X[ID,2,3,Iz,IF] - X[ID,1,3,Iz,IF] 
  end
end  

