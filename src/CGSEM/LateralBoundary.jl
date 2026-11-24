# Open lateral boundaries in x for Cartesian grids (--BoundaryWE="Open"), flow from west to east.
# Inflow (western boundary): the whole state keeps its initial values (zero tendency), and
# after every step it is reset to the stored initial state, which removes the drift caused by
# the implicit vertical solve (OpenBoundaryReset!). A free density at the inflow (one-sided
# strong form derivative) is unstable.
# Outflow (eastern boundary): the grid has one element column more east of the outflow
# x = xOut (ghost elements, added in the driver). Before every right hand side and after every
# step the ghost nodes are filled with the data mirrored at xOut, so the outflow node is an
# interior node with zero x derivative of the whole state after DSS (zero gradient, as at the
# outflow face of a staggered finite volume scheme). The ghost nodes are not part of the domain.
struct OpenBoundary{FT,IT1,AT3}
  xOut::FT
  Col::IT1
  Q0::AT3
  GCol::IT1
  MCol::IT1
end

function OpenBoundary(backend,FE,Metric,U,xOut)
  FT = eltype(U)
  XH = Array(Metric.X)
  GlobH = Array(FE.Glob)
  xMin = minimum(XH[:,1,1,1,:])
  xMax = maximum(XH[:,1,1,1,:])
  tol = 1.0e-6 * (xMax - xMin)
  key(x,y) = (round(Int,x / tol),round(Int,y / tol))
  Col = Int[]
  Phys = Dict{Tuple{Int,Int},Int}()
  for IF in axes(GlobH,2), ID in axes(GlobH,1)
    ic = GlobH[ID,IF]
    xb = XH[ID,1,1,1,IF]
    if abs(xb - xMin) < tol && !(ic in Col)
      push!(Col,ic)
    end
    if xb < xOut + tol
      Phys[key(xb,XH[ID,1,2,1,IF])] = ic
    end
  end
  # ghost node -> physical node mirrored at xOut
  GCol = Int[]
  MCol = Int[]
  for IF in axes(GlobH,2), ID in axes(GlobH,1)
    ic = GlobH[ID,IF]
    xb = XH[ID,1,1,1,IF]
    if xb > xOut + tol && !(ic in GCol)
      push!(GCol,ic)
      push!(MCol,Phys[key(2 * xOut - xb,XH[ID,1,2,1,IF])])
    end
  end
  UH = Array(U)
  Q0H = zeros(eltype(UH),size(UH,2),length(Col),size(UH,4))
  for (IB,ic) in enumerate(Col), iv in axes(UH,4), iz in axes(UH,2)
    Q0H[iz,IB,iv] = UH[1,iz,ic,iv]
  end
  ColD = KernelAbstractions.zeros(backend,Int,length(Col))
  copyto!(ColD,Col)
  Q0 = KernelAbstractions.zeros(backend,FT,size(Q0H))
  copyto!(Q0,Q0H)
  GColD = KernelAbstractions.zeros(backend,Int,length(GCol))
  copyto!(GColD,GCol)
  MColD = KernelAbstractions.zeros(backend,Int,length(MCol))
  copyto!(MColD,MCol)
  return OpenBoundary{FT,typeof(ColD),typeof(Q0)}(FT(xOut),ColD,Q0,GColD,MColD)
end

@kernel inbounds = true function InflowBoundaryKernel!(F,@Const(Col))

  Iz,IB = @index(Global, NTuple)
  NB = @uniform @ndrange()[2]

  if IB <= NB
    ic = Col[IB]
    for iv = 1 : size(F,4)
      F[1,Iz,ic,iv] = eltype(F)(0)
    end
  end
end

@kernel inbounds = true function GhostMirrorKernel!(U,@Const(GCol),@Const(MCol))

  Iz,IB = @index(Global, NTuple)
  NB = @uniform @ndrange()[2]

  if IB <= NB
    ic = GCol[IB]
    jc = MCol[IB]
    for iv = 1 : size(U,4)
      U[1,Iz,ic,iv] = U[1,Iz,jc,iv]
    end
  end
end

function GhostMirror!(U,BC::OpenBoundary,NumberThreadGPU)
  backend = get_backend(U)
  Nz = size(U,2)
  group = (Nz,max(1,div(NumberThreadGPU,Nz)))
  KGhostMirrorKernel! = GhostMirrorKernel!(backend,group)
  KGhostMirrorKernel!(U,BC.GCol,BC.MCol,ndrange=(Nz,length(BC.GCol)))
end
GhostMirror!(U,BC::Nothing,NumberThreadGPU) = nothing

function OpenBoundary!(F,U,p,BC::OpenBoundary,NumberThreadGPU)
  backend = get_backend(F)
  Nz = size(F,2)
  group = (Nz,max(1,div(NumberThreadGPU,Nz)))
  NB = length(BC.Col)
  KInflowBoundaryKernel! = InflowBoundaryKernel!(backend,group)
  KInflowBoundaryKernel!(F,BC.Col,ndrange=(Nz,NB))
end

@kernel inbounds = true function InflowResetKernel!(U,@Const(Col),@Const(Q0))

  Iz,IB = @index(Global, NTuple)
  NB = @uniform @ndrange()[2]

  if IB <= NB
    ic = Col[IB]
    for iv = 1 : size(U,4)
      U[1,Iz,ic,iv] = Q0[Iz,IB,iv]
    end
  end
end

function OpenBoundaryReset!(U,BC::OpenBoundary,NumberThreadGPU)
  backend = get_backend(U)
  Nz = size(U,2)
  NB = length(BC.Col)
  group = (Nz,max(1,div(NumberThreadGPU,Nz)))
  KInflowResetKernel! = InflowResetKernel!(backend,group)
  KInflowResetKernel!(U,BC.Col,BC.Q0,ndrange=(Nz,NB))
  GhostMirror!(U,BC,NumberThreadGPU)
end
OpenBoundaryReset!(U,BC::Nothing,NumberThreadGPU) = nothing
