export JuliaDevice="GPU"
export JuliaGPU="CUDA"
export UCX_ERROR_SIGNALS=""
srun -n 1 gpu_wrapper.sh -n 1 -e "julia --project TestKernels/testKernels.jl"
