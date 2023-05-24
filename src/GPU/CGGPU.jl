  using Metal

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

