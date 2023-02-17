function Plot2DC(c,Fe,xP,zP,Name)
  Nx = size(c,1)
  Nz = size(c,2)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
  cPlot = zeros(Nx*OrdPolyX,Nz*OrdPolyZ)

  for ix = 1 : Nx
    for iz = 1 : Nz  
#     @show c[ix,iz,:,:]  
      cPlot[(ix-1)*OrdPolyX + 1 : ix*OrdPolyX, (iz-1)*OrdPolyZ + 1 : iz*OrdPolyZ] =
        Fe.IntXF2cE*c[ix,iz,:,:]*Fe.IntZC2cE' 
#     @show cPlot[(ix-1)*OrdPolyX + 1 : ix*OrdPolyX, (iz-1)*OrdPolyZ + 1 : iz*OrdPolyZ]  
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
  fig = Figure()
  ax = Axis(fig[1, 1])
  wireframe!(ax, gb_mesh)
  fig
  save(Name*"11.png", fig)
  fig = Figure()
  ax = Axis(fig[1, 1])
  hm = heatmap!(ax, cPlot)
  Colorbar(fig[1, 2], hm)
  fig
  save(Name*".png", fig)
end
