function unstructured_vtkSphere(U,CG,Global, filename::String)
  NF = Global.Grid.NumFaces
  nz = Global.Grid.nz
  OrdPrint = Global.Output.OrdPrint

  Npts = 8 * NF * nz * OrdPrint * OrdPrint
  pts = Array{Float64,2}(undef,3,Npts)
  ipts = 1
  X = zeros(8,3)
  lam=zeros(8,1)
  theta=zeros(8,1)
  z=zeros(8,1)

  for iF = 1 : NumFaces
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
    celltype = VTKCellTypes.VTK_HEX




    ConnectivityList=reshape(1:1:8*NF*OrdPrint*OrdPrint*nz,
    8,NF*OrdPrint*OrdPrint*nz)
    cells = MeshCell[]
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    nf=m.topology.size[3];
    for k in 1 : NF * nz * OrdPrint * OrdPrint
        inds = 1 : 8 + 8 * (k -1)
        push!(cells, MeshCell(celltype, inds))
    end

#   fComp=Array{Array{Float64},1}(undef, length(comp))
#   for l in 1:length(comp)
#       fComp[l]=getElementProperties(p.femType[comp[l]][1],m.meshType,m.geometry.dim,mx,my);
#   end
#   J=Array{Float64,2}(undef,m.geometry.dim,m.topology.dim);
#   dJ=0.0;
#   coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

#   vtk_filename_noext = pwd()*"/output/VTK/"*filename;
#   vtk = vtk_grid(vtk_filename_noext, pts, cells,compress=3)

#   sol=p.solution[tend];
#   for l in 1:length(comp)
#       solc=getfield(sol,comp[l]);
#       if isa(p.degFBoundary[p.femType[comp[l]][1]],degF{1})
#           cvtk=zeros(Float64, nf)
#           for k in 1:nf
#               cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
#               cvtk[k]=dot(fComp[l],cLoc);
#           end
#           vtk_cell_data(vtk, cvtk, name[l])
#       else
#           cvtk=zeros(Float64, m.geometry.dim, nf)
#           for k in 1:nf
#               cLoc=solc[l2g(p.degFBoundary[p.femType[comp[l]][1]], k)]
#               dJ=jacobi!(J,m,k,mx,my,coord);
#               fLoc=(1/dJ)*J*fComp[l]
#               if printSpherical
#                   xyz=transformation(m,coord,mx,my)
#                   lon,lat,r=cart2sphere(xyz[1],xyz[2],xyz[3]);
#                   cvtk[:,k]=velSp(fLoc*cLoc,lon,lat)
#               else
#                   cvtk[:,k]=fLoc*cLoc;
#               end
#           end
#           if size(cvtk,1)==2
#               vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
#               vtk_cell_data(vtk, cvtk[2,:], name[l]*" z")
#               vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],zeros(Float64,nf)), name[l])
#           else
#               vtk_cell_data(vtk, cvtk[1,:], name[l]*" x")
#               vtk_cell_data(vtk, cvtk[2,:], name[l]*" y")
#               vtk_cell_data(vtk, cvtk[3,:], name[l]*" z")
#               vtk_cell_data(vtk, (cvtk[1,:],cvtk[2,:],cvtk[3,:]), name[l])
#           end
#       end
#   end

#   outfiles=vtk_save(vtk);
#   return outfiles::Vector{String}
end
