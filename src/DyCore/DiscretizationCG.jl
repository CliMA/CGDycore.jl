function DiscretizationCG(backend,FT,Jacobi,CG,Exchange,Global,zs) 
# Discretization
  Grid = Global.Grid
  OP=CG.OrdPoly+1
  OPZ=CG.OrdPolyZ+1
  nz=Grid.nz
  NF=Grid.NumFaces

  nQuad = OP * OP
  Metric = MetricStruct{FT}(backend,nQuad,OPZ,Global.Grid.NumFaces,nz)
  dXdxI = zeros(3,3,OPZ,nQuad,nz,NF)
  nS = zeros(nQuad,3,NF)
  FS = zeros(nQuad,NF)
  J = zeros(nQuad,OPZ,nz,NF)
  X = zeros(nQuad,OPZ,3,nz,NF)

  for iF = 1 : NF
    for iz = 1 : nz
      zI = [Grid.z[iz],Grid.z[iz+1]]
      @views (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz) = Jacobi(CG,Grid.Faces[iF],zI,Grids.Topo,Grid.Topography,zs[:,:,iF])
      @views @. X[:,:,:,iz,iF] = X_Fz
      @views @. J[:,:,iz,iF] = J_Fz
      @views @. dXdxI[:,:,:,:,iz,iF] = dXdxI_Fz
      if iz == 1
        #   Surface normal
        @views @. FS[:,iF] = sqrt(dXdxI_Fz[3,1,1,:] * dXdxI_Fz[3,1,1,:] +
          dXdxI_Fz[3,2,1,:] * dXdxI_Fz[3,2,1,:] + dXdxI_Fz[3,3,1,:] * dXdxI_Fz[3,3,1,:])
        @views @. nS[:,1,iF] = dXdxI_Fz[3,1,1,:] / FS[:,iF]
        @views @. nS[:,2,iF] = dXdxI_Fz[3,2,1,:] / FS[:,iF]
        @views @. nS[:,3,iF] = dXdxI_Fz[3,3,1,:] / FS[:,iF]
      end
    end
  end
  copyto!(Metric.J,J)
  copyto!(Metric.X,X)
  copyto!(Metric.dXdxI,dXdxI)
  copyto!(Metric.nS,nS)
  copyto!(Metric.FS,FS)

  MassCGGPU!(backend,FT,CG,Metric.J,CG.Glob,Exchange,Global)

  Metric.dz = KernelAbstractions.zeros(backend,FT,nz,CG.NumG)
  Metric.zP = KernelAbstractions.zeros(backend,FT,nz,CG.NumG)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU  
  NzG = min(div(NumberThreadGPU,OP*OP),nz)  
  group = (OP*OP,NzG,1)  
  ndrange = (OP*OP,nz,NF)
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

  return (CG,Metric)
end

function DiscretizationCG(backend,FT,Jacobi,CG,Exchange,Global)
  DiscretizationCG(backend,FT,Jacobi,CG,Exchange,Global,
  zeros(CG.OrdPoly+1,CG.OrdPoly+1,Global.Grid.NumFaces))
end  

@kernel function GridSizeSphereKernel!(lat,zP,dz,@Const(X),@Const(Glob),Rad)

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]


  if Iz <= Nz && IF <= NF
    @inbounds ind = Glob[ID,IF]
    @inbounds x = eltype(X)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    @inbounds y = eltype(X)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    @inbounds z = eltype(X)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    r = sqrt(x * x + y * y + z * z)
    if Iz == 1
      @inbounds lat[ind] = asin(z / r)
    end  
    @inbounds zP[Iz,ind] = max(r-Rad, eltype(X)(0))
    @inbounds x = X[ID,1,1,Iz,IF]
    @inbounds y = X[ID,1,2,Iz,IF]
    @inbounds z = X[ID,1,3,Iz,IF]
    r1 = sqrt(x * x + y * y + z * z)
    @inbounds x = X[ID,2,1,Iz,IF]
    @inbounds y = X[ID,2,2,Iz,IF]
    @inbounds z = X[ID,2,3,Iz,IF]
    r2 = sqrt(x * x + y * y + z * z)
    @inbounds dz[Iz,ind] =  r2 - r1
  end
end

@kernel function GridSizeCartKernel!(zP,dz,@Const(X),@Const(Glob))

  ID,Iz,IF = @index(Global, NTuple)

  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]
  
  if Iz <= Nz && IF <= NF
    ind = Glob[ID,IF]  
    @atomic zP[Iz,ind] = eltype(X)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF]) 
    @atomic dz[Iz,ind] = X[ID,2,3,Iz,IF] - X[ID,1,3,Iz,IF] 
  end
end  

