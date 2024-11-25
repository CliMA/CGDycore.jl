#!/bin/bash
#PBS -A UCIT0011
#PBS -N gpu_BaroWave_32Elem
#PBS -q main
#PBS -m n
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=64:mpiprocs=4:ngpus=4:mem=400GB
#PBS -o BaroWave_32Elem

# Use scratch for temporary files to avoid space limits in /tmp

# If you are using zsh as default shell


export MPICH_GPU_SUPPORT_ENABLED=1
export LD_PRELOAD=/opt/cray/pe/mpich/8.1.29/gtl/lib/libmpi_gtl_cuda.so.0

export JuliaDevice="GPU"
export JuliaGPU="CUDA"
export machine="derecho"

mpiexec -n 2 -ppn 2 set_gpu_rank ./Jobs/NHSphere/BaroWaveDrySphere_32Elem
