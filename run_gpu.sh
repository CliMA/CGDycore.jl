#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --account=project_465000863
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=2
#SBATCH --mem=16G
#SBATCH --partition=dev-g
# dev-g max 56 cpus per task

#module load CrayEnv rocm

cat << EOF > select_gpu
#!/bin/bash

export MPICH_GPU_SUPPORT_ENABLED=1
export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

echo '####'

#CPU_BIND="map_cpu:49,57"
srun --cpu-bind=cores ./select_gpu julia --project Examples/testNHSphere.jl \
  --Problem="HeldSuarezMoistSphere" \
  --Device="GPU" \
  --GPUType="AMD" \
  --NumberThreadGPU=512 \
  --FloatTypeBackend="Float32" \
  --NumV=5 \
  --NumTr=2 \
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
  --PrintHours=0 \
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
  --HyperDDivW=2.e15

