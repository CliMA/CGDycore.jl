using CGDycore
using LinearAlgebra
n=5
u=zeros(n)
d=zeros(n)
l=zeros(n)
x=zeros(n)
y=zeros(n)
for i=1:n
  u[i] = 2*i  
  l[i] = i  
  d[i] = 3*i  
  x[i] = i*i
end  
u .= u .- 0
l .= l .+ 0
d .= d .+ 10

ev = l[1:n-1]
dv = d
biLL = Bidiagonal(dv, ev, :L)

ev = u[2:n]
dv = d
biUU = Bidiagonal(dv, ev, :U)

triSS = biUU * biLL

tri = zeros(3,n)
biU = zeros(2,n)
biL = zeros(2,n)
biL[1,:] = d
biL[2,1:n-1] = l[1:n-1]
biU[1,2:n] = u[2:n]
biU[2,:] = d

CGDycore.mulUL!(tri,biU,biL)

y .= triSS*x
@show tri
@show y
@show x

CGDycore.triSolve!(x,tri,y)

@show x
