#!/bin/bash
#SBATCH --job-name=my_gpu_job      # Specify job name
#SBATCH --partition=gpu            # Specify partition name
#SBATCH --gpus=8
#SBATCH --exclusive
#SBATCH --mem=0                    # Request all memory available on all nodes
#SBATCH --time=00:30:00            # Set a limit on the total run time
#SBATCH --mail-type=FAIL           # Notify user by email in case of job failure
#SBATCH --account=bb1143           # Charge resources on this project account
#SBATCH --output=OutRace           # File name for standard output

set -e
ulimit -s 204800

# Check GPUs available for the job
# nvidia-smi

# Check GPUs visible for each task
# srun -l nvidia-smi

export JuliaDevice="GPU"
export JuliaGPU="CUDA"
export UCX_ERROR_SIGNALS=""
srun -n 8 gpu_wrapper.sh -n 8 -e "./Jobs/NHSphere/BaroWaveDrySphere_128Elem"
