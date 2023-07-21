
function DivRhoGradE!(F,cC,RhoC,Fe,dXdxI,J,ThreadCache)
  @unpack TCacheC1, TCacheC2, TCacheC3, TCacheC4 = ThreadCache
  Nz = size(F,3)
  N = size(F,1)
  D = Fe.DS
  DW = Fe.DW

  tempx = TCacheC1[Threads.threadid()]
  tempy = TCacheC2[Threads.threadid()]

  @. F = 0
  @inbounds for iz = 1 : Nz
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N  
        Dxc = 0
        Dyc = 0
        @inbounds for k = 1 : N  
          Dxc += D[i,k] * (cC[k,j,iz] / RhoC[k,j,iz])
          Dyc += D[j,k] * (cC[i,k,iz] / RhoC[i,k,iz])
        end
        GradDx = ((dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * Dxc + 
          (dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * Dyc) / (J[i,j,1,iz] + J[i,j,2,iz])
        GradDy = ((dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * Dxc + 
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * Dyc) / (J[i,j,1,iz] + J[i,j,2,iz])
        tempx[i,j] = (dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * GradDx + 
          (dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * GradDy
        tempy[i,j] = (dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * GradDx + 
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * GradDy
      end
    end  
    @inbounds for j = 1 : N
      @inbounds for i = 1 : N  
        @inbounds for k = 1 : N  
          F[i,j,iz] += DW[i,k] * tempx[k,j]
          F[i,j,iz] += DW[j,k] * tempy[i,k] 
        end
      end
    end  
  end    
end    


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
  @tullio FRhoTrC[ix,iy] += D[ix,j] * U[j,iy]
  @tullio FRhoTrC[ix,iy] += D[i,iy] * U[i,ix]
end 
