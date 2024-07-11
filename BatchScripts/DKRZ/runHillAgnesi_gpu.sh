export JuliaDevice="GPU"
export JuliaGPU="CUDA"
export UCX_ERROR_SIGNALS=""
Job=`cat Jobs/NHCart/JobNHHillAgnesiCart`
srun -n 1 gpu_wrapper.sh -n 1 -e "Job"
