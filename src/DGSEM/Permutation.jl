function PermutationO(M,nz)
#Permutation
  N = M * nz 
  p = zeros(Int,3*N)
  ii = 0
  @inbounds for iz = 1 : nz 
    @inbounds for iv = 1 : 1
      @inbounds for k = 1 : M  
        ii += 1 
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end  
    end  
    @inbounds for iv = [2 3] 
#     @inbounds for k = 2 : M - 1
      @inbounds for k = 1 : M - 1
        ii += 1 
        p[ii] = k + (iz - 1) * M + (iv - 1) * N
      end  
    end  
  end  
  ivw = 2
  ivTh = 3
  @inbounds for iz = 1 : nz 
#   ii += 1 
#   p[ii] = 1 + (iz-1) * M  + (ivw - 1) * N
#   ii += 1 
#   p[ii] = 1 + (iz-1) * M  + (ivTh - 1) * N
    ii += 1 
    p[ii] = M + (iz - 1) * M + (ivTh - 1) * N
    ii += 1 
    p[ii] = M + (iz - 1) * M + (ivw - 1) * N
  end  
  return p
end  


function permutation_jacobian_andres(nnodes, nelems)
nvars = 3      
idx = 0
vertical_nodes = nnodes * nelems
perm = zeros(Int,nvars*vertical_nodes)
for elem in 1:nelems
    for node in 1:nnodes
        for v in 1:nvars
            idx += 1
            perm[idx] = Int((nnodes * nelems) * (v - 1) + nnodes * (elem - 1) + node)
        end
    end
end
return perm
end

function perm_to_condensed(nnodes, nelems)
nvars = 3
vertical_nodes = nnodes * nelems
perm = zeros(Int,nvars*vertical_nodes)
size_block_diagonal = (nnodes-1)*nvars*(nelems-1) + nnodes*nvars
idx = 0
idx_block_diagonal = 0
idx_condensed = 0
for elem in 1:nelems-1
    for node in 1:nnodes
        for v in 1:nvars
            idx += 1
            if node == nnodes
            idx_condensed += 1
            perm[idx_condensed+size_block_diagonal] = idx
            else
            idx_block_diagonal+=1
            perm[idx_block_diagonal] = idx
            end
        end
    end
end
## Last element
for node in 1:nnodes
    for v in 1:nvars
        idx += 1
        idx_block_diagonal+=1
        perm[idx_block_diagonal] = idx
    end
end
return perm
end

function permutate_variables_locally(nnodes, nelems)
    polydeg = nnodes - 1
    nvars = 3
    perm_vars = zeros(Int, (polydeg+1)*nvars*nelems)
    idx = 0
    #Permutation for first nelems-1 elements in matrix A
    for elem in 1:nelems-1
        idx_A = (elem-1)*nvars*polydeg
        idx_D = (elem-1)*nvars*polydeg + polydeg*2
        for node in 1:polydeg
            for v in 1:nvars
                idx += 1
                if v == 1 || v == 3
                    idx_A += 1
                    perm_vars[idx_A] = idx
                else
                    idx_D += 1
                    perm_vars[idx_D] = idx
                end
            end
        end
    end
    idx_A = (nelems-1)*nvars*polydeg
    idx_D = (nelems-1)*nvars*polydeg + (polydeg+1)*2
    # One last pass for last element in atrix A
    for node in 1:polydeg+1
        for v in 1:nvars
            idx += 1
            if v == 1 || v == 3
                idx_A += 1
                perm_vars[idx_A] = idx
            else
                idx_D += 1
                perm_vars[idx_D] = idx
            end
        end
    end
    # Just append permutations for matrix D
        for i in idx+1:(polydeg+1)*nvars*nelems
        idx_D +=1
        perm_vars[idx_D] = i
    end
    return perm_vars
end

