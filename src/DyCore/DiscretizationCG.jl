function DiscretizationCG(backend,FT,Jacobi,CG,Global,zs) 
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
  dz = zeros(nz,CG.NumG)
  zP = zeros(nz,CG.NumG)
  J = zeros(nQuad,OPZ,nz,NF)
  X = zeros(nQuad,OPZ,3,nz,NF)

  for iF = 1 : NF
    for iz = 1 : nz
      zI = [Grid.z[iz],Grid.z[iz+1]]
      @views (X_Fz,J_Fz,dXdx_Fz,dXdxI_Fz) = Jacobi(CG,Grid.Faces[iF],zI,Topo,Grid.Topography,zs[:,:,iF])
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

  (M,MW,MMass) = MassCG(CG,J,CG.Glob,Global)
  Global.latN = zeros(CG.NumG)
  latN = Global.latN
  @inbounds for iF = 1 : NF
    @inbounds for iQ = 1 : nQuad
      ind = CG.Glob[iQ,iF]
      @inbounds for iz=1:nz
        dz[iz,ind] +=  2.0*(J[iQ,1,iz,iF] + J[iQ,2,iz,iF])^2 / 
        (dXdxI[3,3,1,iQ,iz,iF] + dXdxI[3,3,2,iQ,iz,iF])  / M[iz,ind]
      end
      @inbounds for iz=1:nz
        if Global.Grid.Form == "Sphere"
          r = norm(0.5 .* (X[iQ,1,:,iz,iF] .+ X[iQ,2,:,iz,iF]))
          zP[iz,ind] += max(r-Global.Grid.Rad, 0.0) * 
            (J[iQ,1,1,iF] + J[iQ,2,1,iF]) / M[iz,ind]
        else
          zP[iz,ind] += 0.5 * (X[iQ,1,3,iz,iF] + X[iQ,2,3,iz,iF]) * 
            (J[iQ,1,1,iF] + J[iQ,2,1,iF]) / M[iz,ind]
        end
      end
    end
  end
  ExchangeData!(dz,Global.Exchange)
  ExchangeData!(zP,Global.Exchange)

  copyto!(Metric.dXdxI,dXdxI)
  copyto!(Metric.nS,nS)
  copyto!(Metric.FS,FS)
  Metric.dz = KernelAbstractions.zeros(backend,FT,size(dz))
  copyto!(Metric.dz,dz)
  Metric.zP = KernelAbstractions.zeros(backend,FT,size(zP))
  copyto!(Metric.zP,zP)
  copyto!(Metric.J,J)
  copyto!(Metric.X,X)
  CG.M = KernelAbstractions.zeros(backend,FT,size(M))
  copyto!(CG.M,M)
  CG.MMass = KernelAbstractions.zeros(backend,FT,size(MMass))
  copyto!(CG.MMass,MMass)
  CG.MW = KernelAbstractions.zeros(backend,FT,size(MW))
  copyto!(CG.MW,MW)
  return (CG,Metric)
end

function DiscretizationCG(backend,FT,Jacobi,CG,Global)
  DiscretizationCG(backend,FT,Jacobi,CG,Global,
  zeros(CG.OrdPoly+1,CG.OrdPoly+1,Global.Grid.NumFaces))
end  
