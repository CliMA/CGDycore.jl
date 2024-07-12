#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1 # number of MPI ranks per node
#SBATCH --gres=gpu:1    # GPUs per node; should equal tasks-per-node
#SBATCH --time=8:00:00
#SBATCH --mem-per-gpu=32GB

module purge
module load climacommon

export JuliaDevice="GPU"
export JuliaGPU="AMD"

srun ./Jobs/NHCart/JobNHHillAgnesiCart

