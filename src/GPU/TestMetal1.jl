using Metal

function memset_kernel(array, value)
  i = thread_position_in_grid_1d()
  if i <= length(array)
    @inbounds array[i] = value
  end
  return
end

a = MtlArray{Float32}(undef, 512)
@metal threads=512 grid=2 memset_kernel(a, 42)

# verify
@assert all(isequal(42), Array(a))
