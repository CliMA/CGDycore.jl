#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --account=project_465000863
#SBATCH --time=00:30:00
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
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
srun --cpu-bind=cores ./select_gpu ./Jobs/NHSphere/JobNHHeldSuarezMoistMOSTSphere
