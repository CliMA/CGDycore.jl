function TopographySmoothing!(Height,CG,Exchange,Global)
 
  backend = get_backend(Height)
  FT = eltype(Height)

  Grid = Global.Grid
  DoF = CG.DoF
  N = CG.OrdPoly + 1
  NF = Global.Grid.NumFaces
  X = KernelAbstractions.zeros(backend,FT,DoF,3,NF)
  dXdxI = KernelAbstractions.zeros(backend,FT,2,2,DoF,NF)
  J = KernelAbstractions.zeros(backend,FT,DoF,NF)
  Rad = Grid.Rad
  F = zeros(4,3,NF)
  FGPU = KernelAbstractions.zeros(backend,FT,4,3,NF)
  for iF = 1 : NF
    F[1,1,iF] = Grid.Faces[iF].P[1].x
    F[1,2,iF] = Grid.Faces[iF].P[1].y
    F[1,3,iF] = Grid.Faces[iF].P[1].z
    F[2,1,iF] = Grid.Faces[iF].P[2].x
    F[2,2,iF] = Grid.Faces[iF].P[2].y
    F[2,3,iF] = Grid.Faces[iF].P[2].z
    F[3,1,iF] = Grid.Faces[iF].P[3].x
    F[3,2,iF] = Grid.Faces[iF].P[3].y
    F[3,3,iF] = Grid.Faces[iF].P[3].z
    F[4,1,iF] = Grid.Faces[iF].P[4].x
    F[4,2,iF] = Grid.Faces[iF].P[4].y
    F[4,3,iF] = Grid.Faces[iF].P[4].z
  end
  copyto!(FGPU,F)
  JacobiSphere2GPU!(X,dXdxI,J,CG,F,Rad)
  M = MassCGGPU2(CG,J,Exchange,Global)

  NFG = min(div(512,N*N),NF)
  group = (N, N, NFG)
  ndrange = (N, N, NF)
  KHyperViscHeightKernel! = HyperViscHeightKernel!(backend,group)
  SmoothType = "Hyper"
  FHeight = similar(Height)
  if SmoothType == "Diff"
     SmoothFac=1.0e9
    @inbounds for i=1:10
      @. FHeight = 0
      KHyperViscHeightKernel!(FHeight,Height,CG.DS,CG.DW,dXdxI,J,M,CG.Glob,ndrange=ndrange) 
      KernelAbstractions.synchronize(backend)
      Parallels.ExchangeData!(FHeight,Exchange)
      @. Height += 0.5 * SmoothFac * FHeight
      @. FHeight = 0
      KHyperViscHeightKernel!(FHeight,Height,CG.DS,CG.DW,dXdxI,J,M,CG.Glob,ndrange=ndrange) 
      KernelAbstractions.synchronize(backend)
      Parallels.ExchangeData!(FHeight,Exchange)
      @. Height += SmoothFac * FHeight
      @. Height = max(Height,0)
    end
  elseif SmoothType == "Hyper"
    SmoothFac=5.e17
    FHeight1 = similar(Height)
    @inbounds for i=1:10
      @. FHeight1 = 0
      KHyperViscHeightKernel!(FHeight1,Height,CG.DS,CG.DW,dXdxI,J,M,CG.Glob,ndrange=ndrange) 
      KernelAbstractions.synchronize(backend)
      Parallels.ExchangeData!(FHeight1,Exchange)
      @. FHeight = 0
      KHyperViscHeightKernel!(FHeight,FHeight1,CG.DS,CG.DW,dXdxI,J,M,CG.Glob,ndrange=ndrange) 
      KernelAbstractions.synchronize(backend)
      Parallels.ExchangeData!(FHeight,Exchange)
      @. Height -= 0.5 * SmoothFac * FHeight

      @. FHeight1 = 0
      KHyperViscHeightKernel!(FHeight1,Height,CG.DS,CG.DW,dXdxI,J,M,CG.Glob,ndrange=ndrange) 
      KernelAbstractions.synchronize(backend)
      Parallels.ExchangeData!(FHeight1,Exchange)
      @. FHeight = 0
      KHyperViscHeightKernel!(FHeight,FHeight1,CG.DS,CG.DW,dXdxI,J,M,CG.Glob,ndrange=ndrange) 
      KernelAbstractions.synchronize(backend)
      Parallels.ExchangeData!(FHeight,Exchange)
      @. Height -= SmoothFac * FHeight
      @show minimum(Height),maximum(Height)
      @. Height = max(Height,0)
    end
  end    
end


@kernel function HyperViscHeightKernel!(FHeight,@Const(Height),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  I, J, iF   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  FaceTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  NF = @uniform @ndrange()[3]

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  HeightCol = @localmem eltype(FHeight) (N,N, FaceTilesDim)
  HeightCxCol = @localmem eltype(FHeight) (N,N, FaceTilesDim)
  HeightCyCol = @localmem eltype(FHeight) (N,N, FaceTilesDim)
  if IF <= NF
    @inbounds HeightCol[I,J,iF] = Height[ind] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if IF <= NF
    @inbounds Dxc = D[I,1] * HeightCol[1,J,iF]
    @inbounds Dyc = D[J,1] * HeightCol[I,1,iF]
    for k = 2 : N
      @inbounds Dxc += D[I,k] * HeightCol[k,J,iF]
      @inbounds Dyc += D[J,k] * HeightCol[I,k,iF] 
    end 
    @inbounds GradDx = (dXdxI[1,1,ID,IF] * Dxc +
      dXdxI[2,1,ID,IF] * Dyc) / JJ[ID,IF]
    @inbounds GradDy = (dXdxI[1,2,ID,IF] * Dxc +
      dXdxI[2,2,ID,IF] * Dyc) / JJ[ID,IF]
    @inbounds tempx = dXdxI[1,1,ID,IF] * GradDx + dXdxI[1,2,ID,IF] * GradDy
    @inbounds tempy = dXdxI[2,1,ID,IF] * GradDx + dXdxI[2,2,ID,IF] * GradDy
    @inbounds HeightCxCol[I,J,iF] = tempx
    @inbounds HeightCyCol[I,J,iF] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if IF <= NF
    @inbounds DivHeight = DW[I,1] * HeightCxCol[1,J,iF] + DW[J,1] * HeightCyCol[I,1,iF]
    for k = 2 : N
      @inbounds DivHeight += DW[I,k] * HeightCxCol[k,J,iF] + DW[J,k] * HeightCyCol[I,k,iF]
    end
    @inbounds @atomic FHeight[ind] += DivHeight / M[ind]
  end
end

function MassCGGPU2(CG,J,Exchange,Global)
  backend = get_backend(J)
  FT = eltype(J)
  N = CG.OrdPoly + 1
  DoF = CG.DoF
  NF = size(CG.Glob,2)

  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU

  NFG = min(div(NumberThreadGPU,N*N),NF)
  group = (N, N, NFG)
  ndrange = (N, N, NF)

  KMassCG2Kernel! = MassCG2Kernel!(backend,group)
  M = KernelAbstractions.zeros(backend,FT,CG.NumG)
  KMassCG2Kernel!(M,J,CG.Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  Parallels.ExchangeData!(M,Exchange)
  return M
end

@kernel function MassCG2Kernel!(M,@Const(JJ),@Const(Glob))
  I,J,IF = @index(Global, NTuple)

  N = @uniform @groupsize()[1]
  NF = @uniform @ndrange()[3]

  if IF <= NF
    ID = I + (J - 1) * N
    @inbounds ind = Glob[ID,IF]
    @inbounds @atomic M[ind] += JJ[ID,IF] 
  end
end
