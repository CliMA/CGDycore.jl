import CGDycore:
  Thermodynamics, Examples, Parallels, Models, Grids, Surfaces,  Outputs, Integration, FiniteElements, DGVertical, GPU, DyCore, DGSEM
using MPI
using Base
using CUDA
using AMDGPU
using Metal
using KernelAbstractions
using StaticArrays
using ArgParse
using LinearAlgebra
using SparseArrays
using BandedMatrices
using FillArrays

function Recover(JacS,M,nz,Phys,fac)
  M2 = M - 2
  m11 = M + 2 * (M - 2)
  A31 = [spzeros(M2,1) Phys.Grav*sparse(I,M2,M2) spzeros(M2,1) ]
  JacSP11 = spzeros(0,0)
  JacSP12 = spzeros(0,0)
  JacSP21 = spzeros(0,nz*m11)
  for iz = 1 : nz
    A13 = JacS.A13[:,:,iz,1]
    A23 = JacS.A23[:,:,iz,1]
    A32 = JacS.A32[:,:,iz,1]
    AS11=[fac*sparse(I,M,M)  spzeros(M,M2) A13
       spzeros(M2,M) fac*sparse(I,M2,M2) A23
       A31 A32 fac*sparse(I,M2,M2)]
    A31F = collect(A31)   
    SA = I(M2)*fac - (1.0 / fac) * A31F * A13 - (1.0 / fac) *A32 * A23
    DGSEM.LUFull!(SA)
    @show sum(abs.(SA - JacS.SA[:,:,iz,1]))
    JacSP11 = [JacSP11 spzeros(m11*(iz-1),m11)
            spzeros(m11,m11*(iz-1)) AS11]
    B1_1 = JacS.B1_1[iz,1]
    B1_23 = JacS.B1_23[:,:,iz,1]
    B1_4 = JacS.B1_4[iz,1]
    B1m_34 = JacS.B1m_34[:,iz,1]'
    B1p_12 = JacS.B1p_12[:,iz,1]'
    B2_23 = JacS.B2_23[:,:,iz,1]
    B3_14  = JacS.B3_14[:,:,iz,1]
    BS12_l = [spzeros(m11,2) [B1m_34
                            spzeros(m11-1,2)]]
    BS12_r = [[spzeros(M-1,2)
               B1p_12
               spzeros(2*M2,2)] spzeros(m11,2)]
    BS12_c1 = [B1_1
               spzeros(M-1+M2,1)
               B3_14[:,1]] 
    BS12_c23 = [B1_23
                B2_23
                spzeros(M2,2)]
    BS12_c4 = [spzeros(M-1,1)
               B1_4
               spzeros(M2,1)
               B3_14[:,2]]
    BS12_c = [BS12_c1 BS12_c23 BS12_c4]           
    if iz == 1
      JacSP12 = [BS12_c BS12_r spzeros(m11,(nz-2)*4)]
    elseif iz == nz
      JacSP12 = [JacSP12
                 spzeros(m11,(nz-2)*4) BS12_l BS12_c]
    else             
      JacSP12 = [JacSP12
                 spzeros(m11,(iz-2)*4) BS12_l BS12_c BS12_r spzeros(m11,(nz-1-iz)*4)]
    end             
    C14_3 = JacS.C14_3[:,:,iz,1]
    C23_2 = JacS.C23_2[:,:,iz,1]
    CS11 = [spzeros(1,M)
            [Phys.Grav spzeros(1,M-1)]
            [spzeros(1,M-1) Phys.Grav ]
            spzeros(1,M)]
    CS12 = [spzeros(1,M2)
            C23_2
            spzeros(1,M2)]
    CS13 = [C14_3[1,:]'
            spzeros(2,M2)
            C14_3[2,:]']
    JacSP21 = [JacSP21
               [spzeros(4,(iz-1)*m11) CS11 CS12 CS13 spzeros(4,(nz-iz)*m11)]]

  end
  B22 = BandedMatrix{Float64}(undef, (4*nz, 4*nz), (3, 3))
  @. B22.data = JacS.SchurBand[:,:,1]
  JacSP22 = sparse(B22)
  return JacSP11,JacSP12,JacSP21,JacSP22
end


backend = CPU()
FTB = Float64
Parallel = true

MPI.Init()
comm = MPI.COMM_WORLD
Proc = MPI.Comm_rank(comm) + 1
ProcNumber = MPI.Comm_size(comm)
ParallelCom = DyCore.ParallelComStruct()
ParallelCom.Proc = Proc
ParallelCom.ProcNumber  = ProcNumber

Problem = "HillAgnesiXCart"
Param = Examples.Parameters(FTB,Problem)
# Physical parameters
Phys = DyCore.PhysParameters{FTB}()

# ModelParameters
Model = DyCore.ModelStruct{FTB}()

# Grid
Boundary = Grids.Boundary()
Boundary.WE = "Period"
Boundary.SN = "Period"
Boundary.BT = "FreeSlip"
nx = 2
ny = 2
Lx = 1000.0
Ly = 1000.0
x0 = 0.0
y0 = 0.0
nz = 5
OrdPoly = 4
OrdPolyZ = 4
M = OrdPolyZ + 1
N = M * nz
Grid, CellToProc = Grids.InitGridCart(backend,FTB,OrdPoly,nx,ny,Lx,Ly,x0,y0,Boundary,nz,Model,ParallelCom)

OrdPrint = 4
OrdPrintZ = 4
DG = FiniteElements.DGQuad{FTB}(backend,OrdPoly,OrdPolyZ,OrdPrint,OrdPrintZ,Grid,ParallelCom.Proc)


DG.NumG = 1
DG.NumI = 1
fac = 0.5
U = ones(M,nz,DG.NumG,5)
H = 10000.0
dzLoc = H / nz
dz = ones(nz,DG.NumG) * dzLoc

Profile = Examples.StratifiedExample()(Param,Phys)
time = 0.0
for iz = 1 : nz
  z0 = (iz - 1) * dzLoc  
  z1 = iz * dzLoc  
  for k = 1 : M
    zLoc = 0.5 * ((1.0 - DG.xwZ[k]) * z0 + (1.0 + DG.xwZ[k]) * z1)  
    xS = SVector{3}(0.0,0.0,zLoc)
    RhoP,_,_,_,ThP= Profile(xS,time)
    @views @. U[k,iz,:,1] = RhoP
    @views @. U[k,iz,:,5] = RhoP * ThP
  end
end

dSdS,dSdM,dMdS,dMdM = DGSEM.InitJacDG(DG,nz,Param)
JacS = DGSEM.JacDGVert{FTB}(backend,M,nz,DG.NumI)
JacGPU = DGSEM.JacDGVert{FTB}(backend,M,nz,DG.NumI)

fac = 2.0
DGSEM.FillJacDGVertOld!(JacS,U,DG,dz,fac,Phys,Param)
DGSEM.FillJacDGVert!(JacGPU,U,DG,dz,fac,Phys,Param)
BP11,BP12,BP21,BP22= Recover(JacS,M,nz,Phys,fac)
p = DGSEM.Permutation(M,nz)
pI = invperm(p)
BP = [BP11 BP12
      BP21 BP22]
B = BP[pI,pI]      

@show sum(abs.(JacS.SchurBand-JacGPU.SchurBand))
DGSEM.SchurBoundaryOld!(JacS)
DGSEM.SchurBoundary!(JacGPU)
@show sum(abs.(JacS.SchurBand-JacGPU.SchurBand))
stop

JacVLU, A = DGSEM.JacDG(U,DG,fac,dSdS,dSdM,dMdS,dMdM,dz,Phys)

b = ones(size(U))
@. b[:,:,:,5] *= 2
@. b[:,:,:,4] *= 3
x = similar(b)
@. x = b

bF = zeros(M*nz*3)
for ID = 1 : DG.NumI
  ibF = 0
  @inbounds for iv in [1,4,5]
    @inbounds for iz = 1 : nz
      @inbounds for i = 1 : M
        ibF += 1
        bF[ibF] = x[i,iz,ID,iv]
      end
    end
  end
  ldiv!(JacVLU[ID],bF)
  ibF = 0
  @inbounds for iv in [1,4,5]
    @inbounds for iz = 1 : nz
      @inbounds for i = 1 : M
        ibF += 1
        x[i,iz,ID,iv] = bF[ibF]
      end
    end
  end
end
DGSEM.ldivVertical!(JacS,b)
@show sum(abs.(x-b))

aaa=3
