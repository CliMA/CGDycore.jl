mutable struct Metric
  X::Array{Float64, 7}
  J::Array{Float64, 6}
  JC::Array{Float64, 6}
  dXdx::Array{Float64, 8}
  dXdxI::Array{Float64, 8}
  norm_nSB::Array{Float64, 5}
  nSB::Array{Float64, 6}
  VolSurfB::Array{Float64, 5}
  norm_nST::Array{Float64, 5}
  nST::Array{Float64, 6}
  VolSurfT::Array{Float64, 5}
end

function Metric(Nx,Ny,Nz,xP,yP,zP,Fe)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
# Metric
  X=zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,3)
  J=zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  JC=zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ)
  dXdx=zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,3,3)
  dXdxI=zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,3,3)
  norm_nSB = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1)
  nSB = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,3)
  VolSurfB = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1)
  norm_nST = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1)
  nST = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,3)
  VolSurfT = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1)

  for iz = 1 : Nz
    for ix = 1 : Nx
      for iy = 1 : Ny
        @views (XLoc,JLoc,dXdxLoc,dXdxILoc,norm_nSBLoc,nSBLoc,VolSurfBLoc,
          norm_nSTLoc,nSTLoc,VolSurfTLoc)=JacobiDG3(xP[ix,iy,iz,:,:,:],
          yP[ix,iy,iz,:,:,:],zP[ix,iy,iz,:,:,:],Fe)
        @views @. X[ix,iy,iz,:,:,:,:] = XLoc
        @views @. J[ix,iy,iz,:,:,:] = JLoc
        @views @. dXdx[ix,iy,iz,:,:,:,:,:] = dXdxLoc
        @views @. dXdxI[ix,iy,iz,:,:,:,:,:] = dXdxILoc
        @views @. norm_nSB[ix,iy,iz,:,:] = norm_nSBLoc
        @views @. nSB[ix,iy,iz,:,:,:] = nSBLoc
        @views @. VolSurfB[ix,iy,iz,:,:] = VolSurfBLoc
        @views @. norm_nST[ix,iy,iz,:,:] = norm_nSTLoc
        @views @. nST[ix,iy,iz,:,:,:] = nSTLoc
        @views @. VolSurfT[ix,iy,iz,:,:] = VolSurfTLoc
        for i = 1 : OrdPolyX +1
          for j = 1 : OrdPolyY +1
            @views JC[ix,iy,iz,i,j,:] = Fe.IntZF2C * J[ix,iy,iz,i,j,:]
          end
        end
      end
    end
  end
  return Metric(
    X,
    J,
    JC,
    dXdx,
    dXdxI,
    norm_nSB,
    nSB,
    VolSurfB,
    norm_nST,
    nST,
    VolSurfT,
  )
end  


