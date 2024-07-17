#!/bin/bash
#PBS -A UCIT0011
#PBS -N gpu_aquaplanet_dyamond_helem_256_0M_256gpu_job
#PBS -q main
#PBS -m n
#PBS -l walltime=01:00:00
#PBS -l select=2:ncpus=8:mem=480GB:ngpus=4

# Use scratch for temporary files to avoid space limits in /tmp

# If you are using zsh as default shell
source /glade/u/apps/derecho/23.09/spack/opt/spack/lmod/8.7.24/gcc/7.5.0/c645/lmod/lmod/init/zsh

export MODULEPATH="/glade/campaign/univ/ucit0011/ClimaModules-Derecho:$MODULEPATH"
#module --force purge
module load climacommon/2024_05_27
module load julia-preferences/2024_02_20
module load julia/1.10.4
module load mpiwrapper/2024_05_27

export JuliaDevice="GPU"
export JuliaGPU="CUDA"

mpiexec -n 4 -ppn 4 ./Jobs/NHCart/JobNHHillAgnesiCart
