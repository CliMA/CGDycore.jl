using BenchmarkTools
using LinearAlgebra

function progonka!(u,a,b,c,f)
    N = size(f,1)
    α=similar(f)
    β=similar(f)
    α[2] = -c[1]/b[1]
    β[2] = f[1]/b[1]
    @inbounds @fastmath for i in 2:N-1 #Forward Sweep
        α[i+1]=-c[i]/(a[i-1]*α[i]+b[i])
        β[i+1]=(f[i]-a[i-1]*β[i])/(a[i-1]*α[i]+b[i])
    end
    u[N]=(f[N]-a[N-1]*β[N])/(a[N-1]*α[N]+b[N])
    @inbounds @fastmath for i in N-1:-1:1 #Backward Sweep
        u[i]=α[i+1]*u[i+1]+β[i+1]
    end
    u
end

progonka(a,b,c,f)=progonka!(similar(f),a,b,c,f)


function tprogonka()
    N=1000
    a=-rand(N-1)
    c=-rand(N-1)
    b=2.0.+rand(N)
    f=rand(N)
   
    A=Tridiagonal(a,b,c)
    utri=A\f
    upro=progonka(a,b,c,f)
    @assert utri≈upro
    @btime utri=$A\$f
    @btime upro=progonka($a,$b,$c,$f)
    nothing
end
