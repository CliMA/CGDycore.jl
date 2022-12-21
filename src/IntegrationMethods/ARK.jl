struct ARS343 <: AbstractIMEXARKTableau end
function tableau(::ARS343)
    γ = 0.4358665215084590
    a42 = 0.5529291480359398
    a43 = a42
    b1 = -3/2 * γ^2 + 4 * γ - 1/4
    b2 =  3/2 * γ^2 - 5 * γ + 5/4
    a31 = (1 - 9/2 * γ + 3/2 * γ^2) * a42 +
        (11/4 - 21/2 * γ + 15/4 * γ^2) * a43 - 7/2 + 13 * γ - 9/2 * γ^2
    a32 = (-1 + 9/2 * γ - 3/2 * γ^2) * a42 +
        (-11/4 + 21/2 * γ - 15/4 * γ^2) * a43 + 4 - 25/2 * γ + 9/2 * γ^2
    a41 = 1 - a42 - a43
    return IMEXARKTableau(;
        a_exp = @SArray([
            0   0   0   0;
            γ   0   0   0;
            a31 a32 0   0;
            a41 a42 a43 0;
        ]),
        b_exp = @SArray([0, b1, b2, γ]),
        a_imp = @SArray([
            0 0       0  0;
            0 γ       0  0;
            0 (1-γ)/2 γ  0;
            0 b1      b2 γ;
        ])
    )
end


function step_u!(integrator, cache::NewIMEXARKCache)
    (; u, p, t, dt, sol, alg) = integrator
    (; f) = sol.prob
    (; T_lim!, T_exp!, T_imp!, lim!, dss!, stage_callback!) = f
    (; tab, newtons_method) = alg
    (; a_exp, b_exp, a_imp, b_imp, c_exp, c_imp) = tab
    (; U, T_lim, T_exp, T_imp, temp, γ, newtons_method_cache) = cache
    s = length(b_exp)

    if !isnothing(T_imp!)
        update!(
            newtons_method,
            newtons_method_cache,
            NewTimeStep(t),
            jacobian -> isnothing(γ) ?
                error(
                    "The tableau does not specify a unique value of γ for the \
                     duration of each time step; do not update based on the \
                     NewTimeStep signal when using this tableau."
                ) : T_imp!.Wfact(jacobian, u, p, dt * γ, t),
        )
    end

    for i in 1:s
      t_exp = t + dt * c_exp[i]
      t_imp = t + dt * c_imp[i]

      @. U[i] = u

      for j in 1:(i - 1)
        @. U[i] += dt * a_exp[i, j] * T_exp[j]
      end

      for j in 1:(i - 1)
        @. U[i] += dt * a_imp[i, j] * T_imp[j]
      end


      if a_imp[i, i] >0
        @. temp = U[i]
        for iterN = 1 : 2
           T_imp!(residual, U[i], p, t_imp)
           @. residual = temp + dt * a_imp[i, i] * residual - U[i]
           # Solve U[i] = sss
        end    
     end
     T_exp!(T_exp[i], U[i], p, t_exp)
     if a_imp[i, i] >0
       @. T_imp[i] = (U[i] - temp) / (dt * a_imp[i, i])
     else
       T_imp!(T_imp[i], U[i], p, t_imp)
     end  
    end

    t_final = t + dt

    for j in 1:s
     @. u += dt * b_exp[j] * T_exp[j]
    end
    for j in 1:s
       @. u += dt * b_imp[j] * T_imp[j]
    end
    return u
end
