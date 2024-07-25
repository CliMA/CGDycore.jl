using CUDA, AMDGPU
using KernelAbstractions
using KernelAbstractions: @atomic
using KernelAbstractions.Extras
using ImageShow, ImageIO


function index_fun(arr; backend=get_backend(arr))
	out = similar(arr)
	fill!(out, 0)
	kernel! = my_kernel!(backend)
	kernel!(out, arr, ndrange=(size(arr, 1), size(arr, 2)))
        KernelAbstractions.synchronize(backend)
        @time for i = 1 : 200
	    kernel!(out, arr, ndrange=(size(arr, 1), size(arr, 2)))
            KernelAbstractions.synchronize(backend)
        end    
	kernel1! = my_kernel1!(backend)
	kernel1!(out, arr, ndrange=(size(arr, 1), size(arr, 2)))
        KernelAbstractions.synchronize(backend)
        @time for i = 1 : 200
	    kernel1!(out, arr, ndrange=(size(arr, 1), size(arr, 2)))
            KernelAbstractions.synchronize(backend)
        end    
	kernel2! = my_kernel2!(backend)
	kernel2!(out, arr, ndrange=(size(arr, 1), size(arr, 2)))
        KernelAbstractions.synchronize(backend)
        @time for i = 1 : 200
	    kernel2!(out, arr, ndrange=(size(arr, 1), size(arr, 2)))
            KernelAbstractions.synchronize(backend)
        end    
	return out
end

@kernel function my_kernel!(out, arr)
	i, j = @index(Global, NTuple)
	for k in 1:size(out, 1)
             @atomic :monotonic out[k, i] += arr[i, j]
	end
end

@kernel function my_kernel1!(out, arr)
        i, j = @index(Global, NTuple)
        for k in 1:size(out, 1)
             @atomic out[k, i] += arr[i, j]
        end
end

@kernel function my_kernel2!(out, arr)
        i, j = @index(Global, NTuple)
        for k in 1:size(out, 1)
             out[k, i] += arr[i, j]
        end
end

JuliaDevice = get(ENV, "JuliaDevice", "CPU")
JuliaGPU = get(ENV, "JuliaGPU", "CUDA")

if JuliaDevice == "CPU"
  backend = CPU()
elseif JuliaDevice == "GPU"
  if JuliaGPU == "CUDA"
    backend = CUDABackend()
    CUDA.allowscalar(false)
#   CUDA.device!(MPI.Comm_rank(MPI.COMM_WORLD))
  elseif JuliaGPU == "AMD"
    backend = ROCBackend()
    AMDGPU.allowscalar(false)
  elseif JuliaGPU == "Metal"
    backend = MetalBackend()
    Metal.allowscalar(true)
  end
else
  backend = CPU()
end

FT = Float32

imgGPU = KernelAbstractions.zeros(backend, FT, (500, 500));
img = zeros(FT, (500, 500));
img[100:200, 100:200] .= 1;
img[350:450, 350:450] .= 2;
copyto!(imgGPU,img)

out = Array(index_fun(imgGPU));
