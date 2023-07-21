using CGDycore
using MPI
using Base


function DivRhoGradE1!(F,cC,RhoC,D,DW,dXdxI,J)

  Nz = size(F,3)
  N = size(F,1)

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
        tempx = (dXdxI[i,j,1,iz,1,1] + dXdxI[i,j,2,iz,1,1]) * GradDx +
          (dXdxI[i,j,1,iz,1,2] + dXdxI[i,j,2,iz,1,2]) * GradDy
        tempy = (dXdxI[i,j,1,iz,2,1] + dXdxI[i,j,2,iz,2,1]) * GradDx +
          (dXdxI[i,j,1,iz,2,2] + dXdxI[i,j,2,iz,2,2]) * GradDy
        for k = 1 : N
          F[k,j,iz] += DW[k,i] * tempx
          F[i,k,iz] += DW[k,j] * tempy
        end
      end
    end
  end
end

Nz = 10
OrdPoly = 4

(DW,DS,) = CGDycore.DerivativeMatrixSingle(OrdPoly)

F = zeros(Float64,OrdPoly+1,OrdPoly+1,Nz)
cC = rand(Float64,OrdPoly+1,OrdPoly+1,Nz)
RhoC = ones(Float64,OrdPoly+1,OrdPoly+1,Nz)
dXdxI = rand(OrdPoly+1,OrdPoly+1,2,Nz,3,3)
J = ones(OrdPoly+1,OrdPoly+1,2,Nz)

DivRhoGradE1!(F,cC,RhoC,DS,DW,dXdxI,J)
