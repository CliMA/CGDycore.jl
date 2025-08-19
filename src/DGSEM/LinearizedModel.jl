function FcnLin!(F,U,DG,dx,dz)

  nx = length(dx)
  nz = lenghth(dz)
  @views Th = Aux[:,:,:,:,1]
  @views pTh = Aux[:,:,:,:,2]

  @views Rho = U[:,:,:,:,1]
  @views Rhou = U[:,:,:,:,2]
  @views Rhow = U[:,:,:,:,3]
  @views RhoTh = U[:,:,:,:,4]

  @views FRho = F[:,:,:,:,1]
  @views FRhou = F[:,:,:,:,2]
  @views FRhow = F[:,:,:,:,3]
  @views FRhoTh = F[:,:,:,:,4]

  DW = DG.DW

  #Volume kernel

  for ix = 1 : nx
    for iz = 1 : nz
      FRho[:,:,ix,iz] = (1 / dx[ix]) * DW * Rhou[:,:,ix,iz] + (1 / dz[iz]) * Rhow[:,:,ix,iz] * DW' 
      FRhou[:,:,ix,iz] = (1 / dx[ix]) * DW * (pTh[:,:,ix,iz] .* RhoTh[:,:,ix,iz])
      FRhow[:,:,ix,iz] = (1 / dz[iz]) * (pTh[:,:,ix,iz] .* RhoTh[:,:,ix,iz]) * DW'
      FRho[:,:,ix,iz] = (1 / dx[ix]) * DW * (Rhou[:,:,ix,iz] .* Th[:,:,ix,z]) + 
        (1 / dz[iz]) * (Rhow[:,:,ix,iz] .* Th[:,:,ix,z]) * DW' 
    end
  end   
end
