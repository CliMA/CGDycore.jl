function vadd(a, b, c)
  i = thread_position_in_grid_1d()
  c[i] = a[i] + b[i]
  return
end

function vmult(a, b, c)
  i = thread_position_in_grid_1d()
  c[i] = a[i] * b[i]
  return
end

function test(n=10, k=20, m=30; T=Int16, U=Float32)
    a = rand(T, n, k)
    b = rand(T, k, m)
    c = U.(a) * U.(b)   # CPU can't do mixed-precision (or at least,
                        # not without losing precision)

    d_a = MtlArray(a)
    d_b = MtlArray(b)
    d_c = MtlArray(similar(c))
    mul!(d_c, d_a, d_b)

    display(d_c)
    display(c)

    Array(d_c) ≈ c
end
