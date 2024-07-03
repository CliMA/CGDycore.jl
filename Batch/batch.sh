!/bin/bash
#SBATCH --account=project_465000863
#SBATCH --partition=small-g
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=0

module use /appl/local/csc/modulefiles
module load julia
module load julia-mpi
module load julia-amdgpu
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project Examples/testNHSphere.jl \
  --Problem="HeldSuarezMoistSphere" \
  --Device="CPU" \
  --GPUType="Metal" \
  --NumberThreadGPU=512 \
  --FloatTypeBackend="Float32" \
  --NumV=5 \
  --NumTr=2 \
  --TkePos=0 \
  --RhoVPos=6 \
  --RhoCPos=7 \
  --ProfpBGrd="" \
  --ProfRhoBGrd="" \
  --Source=true \
  --Forcing=true \
  --Curl=false \
  --ModelType="VectorInvariant" \
  --VerticalDiffusion=true \
  --SurfaceFlux=true \
  --SurfaceScheme="" \
  --Coriolis=true \
  --Upwind=true \
  --HorLimit=false \
  --Equation="CompressibleShallow" \
  --State="Moist" \
  --Microphysics=true \
  --TypeMicrophysics="SimpleMicrophysics" \
  --Buoyancy=true \
  --Damping=true \
  --StrideDamp=10000 \
  --Relax=1.0e-3 \
  --Decomp="EqualArea" \
  --SimSeconds=0 \
  --SimDays=10 \
  --PrintSeconds=0 \
  --PrintMinutes=0 \
  --PrintHours=12 \
  --PrintDays=0 \
  --StartAverageDays=100 \
  --Flat=true \
  --dtau=150 \
  --IntMethod="Rosenbrock" \
  --Table="SSP-Knoth" \
  --TopoS="" \
  --Stretch=true \
  --StretchType="Exp" \
  --GridType="CubedSphere" \
  --nz=64 \
  --nPanel=30 \
  --H=45000.0 \
  --OrdPoly=3 \
  --HyperVisc=true \
  --HyperDCurl=2.e15 \
  --HyperDGrad=2.e15 \
  --HyperDDiv=2.e15 \
  --HyperDDivW=0.e15
