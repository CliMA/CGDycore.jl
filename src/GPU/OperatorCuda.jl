
function DivRhoTrCellCUDA!(FRhoTrC,uC,vC,RhoTrC,Fe,dXdxI,U,V)
  OrdPoly = Fe.OrdPoly
  D = Fe.DS

  ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
  iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
  U[ix,iy] = -RhoTrC[ix,iy] * ((dXdxI[ix,iy,1,1,1] + dXdxI[ix,iy,2,1,1]) * uC[ix,iy] +
    (dXdxI[ix,iy,1,1,2] + dXdxI[ix,iy,2,1,2]) * vC[ix,iy])
  V[ix.iy] = -RhoTrC[ix,iy] * ((dXdxI[ix,iy,1,2,1] + dXdxI[ix,iy,2,2,1]) * uC[ix,iy] +
    (dXdxI[ix,iy,1,2,2] + dXdxI[ix,iy,2,2,2]) * vC[ix,iy])
  sync_threads()
  ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
  iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
  @tullio FRhoTrC[ix,iy] += D[i,iy] * U[i,iy]
  @tullio FRhoTrC[ix,iy] += D[ix,j] * U[ix,j]
end 
