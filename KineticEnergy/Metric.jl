mutable struct Metric
  X::Array{Float64, 5}
  J::Array{Float64, 4}
  JC::Array{Float64, 4}
  dXdx::Array{Float64, 6}
  dXdxI::Array{Float64, 6}
  norm_nSB::Array{Float64, 3}
  nSB::Array{Float64, 4}
  VolSurfB::Array{Float64, 3}
  norm_nST::Array{Float64, 3}
  nST::Array{Float64, 4}
  VolSurfT::Array{Float64, 3}
end

function Metric(Nx,Nz,xP,zP,Fe)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyZ = Fe.OrdPolyZ
# Metric
  X=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1,2)
  J=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1)
  JC=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ)
  dXdx=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1,2,2)
  dXdxI=zeros(Nx,Nz,OrdPolyX+1,OrdPolyZ+1,2,2)
  norm_nSB = zeros(Nx,Nz,OrdPolyX+1)
  nSB = zeros(Nx,Nz,OrdPolyX+1,2)
  VolSurfB = zeros(Nx,Nz,OrdPolyX+1)
  norm_nST = zeros(Nx,Nz,OrdPolyX+1)
  nST = zeros(Nx,Nz,OrdPolyX+1,2)
  VolSurfT = zeros(Nx,Nz,OrdPolyX+1)

  for iz = 1 : Nz
    for ix = 1 : Nx
      @views (XLoc,JLoc,dXdxLoc,dXdxILoc,norm_nSBLoc,nSBLoc,VolSurfBLoc,
        norm_nSTLoc,nSTLoc,VolSurfTLoc)=JacobiDG2(xP[ix,iz,:,:],zP[ix,iz,:,:],Fe)
      @views @. X[ix,iz,:,:,:] = XLoc
      @views @. J[ix,iz,:,:] = JLoc
      @views @. dXdx[ix,iz,:,:,:,:] = dXdxLoc
      @views @. dXdxI[ix,iz,:,:,:,:] = dXdxILoc
      @views @. norm_nSB[ix,iz,:] = norm_nSBLoc
      @views @. nSB[ix,iz,:,:] = nSBLoc
      @views @. VolSurfB[ix,iz,:] = VolSurfBLoc
      @views @. norm_nST[ix,iz,:] = norm_nSTLoc
      @views @. nST[ix,iz,:,:] = nSTLoc
      @views @. VolSurfT[ix,iz,:] = VolSurfTLoc
      for i = 1 : OrdPolyX +1
        @views JC[ix,iz,i,:] = Fe.IntZF2C * J[ix,iz,i,:]
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


