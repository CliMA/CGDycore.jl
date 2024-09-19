#!/bin/bash
#SBATCH --job-name=my_gpu_job      # Specify job name
#SBATCH --partition=gpu            # Specify partition name
#SBATCH --gpus-per-node=4
#SBATCH --nodes=16   # -> 8GPUs
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=64
#SBATCH --exclusive
#SBATCH --mem=0                    # Request all memory available on all nodes
#SBATCH --time=00:30:00            # Set a limit on the total run time
#SBATCH --mail-type=FAIL           # Notify user by email in case of job failure
#SBATCH --account=bb1143           # Charge resources on this project account
#SBATCH --output=BaroWave_320Elem           # File name for standard output

set -e
ulimit -s 204800

# Check GPUs available for the job
# nvidia-smi

# Check GPUs visible for each task
# srun -l nvidia-smi

export JuliaDevice="GPU"
export JuliaGPU="CUDA"
export machine="levante"
export UCX_ERROR_SIGNALS=""
srun -n 64 gpu_wrapper.sh -n 64 -e "./Jobs/NHSphere/BaroWaveDrySphere_320Elem"
