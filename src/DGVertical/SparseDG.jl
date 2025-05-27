using LinearAlgebra
using SparseArrays

function Jacobian(U,DG)
  
  FT = eltype(U) 
  Nz = size(U,2)
  M = size(U,1)
  N = Nz * M
  RowInd = Int[]
  ColInd = Int[]
  Val = FT[]
  for iZ = 1 : Nz
    for i = 1 : M  
      for j = 1 : M  
        push!(RowInd,i+(iZ-1)*M)  
        push!(ColInd,j+(iZ-1)*M)  
        push!(DG.DS[i,j])
    end
  end  
  if iZ < Nz
    push!(RowInd,1+iZ*M)  
    push!(ColInd,M+(iZ-1)*M)  
    push!(Val,1.0)
    push!(RowInd,iZ*M)  
    push!(ColInd,1+iZ*M)  
    push!(Val,1.0)
  end  
end
Grad = sparse(RowInd, ColInd, Val)
RowInd = Int[]
ColInd = Int[]
Val = Int[]
for iZ = 1 : Nz
  if iZ < Nz
    push!(RowInd,iZ*M)
    push!(ColInd,M+(iZ-1)*M)
    push!(Val,1.0)
    push!(RowInd,iZ*M)
    push!(ColInd,1+iZ*M)
    push!(Val,1.0)
    push!(RowInd,1+iZ*M)
    push!(ColInd,iZ*M)
    push!(Val,1.0)
    push!(RowInd,1+iZ*M)
    push!(ColInd,1+iZ*M)
    push!(Val,1.0)
  end
end

permu = zeros(Int,N)
iInd = 1
shift = (M-2) * Nz
for iZ = 1 : Nz
  for i = 2 : M - 1
    permu[iInd] = i + (iZ-1)*M
    global iInd += 1
  end  
end  
#for iZ = 1 : Nz
#  for i = 2 : M - 1
#    permu[iInd] = i + (iZ-1)*M + N
#    global iInd += 1
#  end  
#end  
for iZ = 1 : Nz
  permu[iInd] = 1+(iZ-1)*M
  global iInd += 1
  permu[iInd] = M + (iZ-1)*M
  global iInd += 1
end
#for iZ = 1 : Nz
#  permu[iInd] = 1+(iZ-1)*M + N
#  global iInd += 1
#  permu[iInd] = M + (iZ-1)*M + N
#  global iInd += 1
#end

GradDiv = Grad*Grad
GradDivP = GradDiv[permu,permu]
GradDivP11 = GradDivP[1:(M-2)*Nz,1:(M-2)*Nz]
GradDivP12 = GradDivP[1:(M-2)*Nz,(M-2)*Nz+1:end]
GradDivP21 = GradDivP[(M-2)*Nz+1:end,1:(M-2)*Nz]
GradDivP22 = GradDivP[(M-2)*Nz+1:end,(M-2)*Nz+1:end]
S = GradDivP22 - 10*GradDivP21 * GradDivP12
SS = GradDivP21 * GradDivP12

Id = sparse(I,N,N)
B = Id - 0.01 * GradDivP
B = Id - 0.01 * GradDiv
Chol = cholesky(B,perm=1:N)
L = sparse(Chol.L)

D = sparse(RowInd, ColInd, Val, N, N)
Id = sparse(I,2*N,2*N)
DP=D[permu,permu]
GradP = Grad[permu,permu]
Jac=[DP GradP
   GradP DP]
#=   
JacP = Jac[permu,permu]   
gamma = 0.1   
B = Id - gamma * JacP   
q = zeros(Int,2*N)
q .= 1 : 2*N
LUB = lu(B,q=q)
CholB = cholesky(B,perm=q)
=#

