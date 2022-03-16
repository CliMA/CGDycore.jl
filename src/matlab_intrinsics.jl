strcmp(s1, s2) = s1 == s2
strcat(s...) = prod(s)
num2str(n) = "$n"
lower(s) = lowercase(s)
repmat(v::AbstractArray,s...) = repeat(v, s...)
repmat(v::Real,s...) = repeat([v], s...)

function atan2(y, x)
    return atan(y, x)
end

function permute(args...)
    permutedims(args...)
end

# spdiags([reshape(-D,nJ,1) reshape(D,nJ,1)],[-1 0] ,nJ,nJ)
function spdiags(mat, opt::Real, m, n)
    return Diagonal(mat)
end
function spdiags(mat, opt::AbstractArray, m, n)
    @assert length(opt) == 2
    # @show size(mat)
    # @show opt
    # error("Debugging spdiags")
    if opt[1] == -1 && opt[2] == 0
        ev = mat[:,1][2:end] # TODO: verify
        dv = mat[:,2]
        return Bidiagonal(dv, ev, :L)
    elseif opt[1] == 0 && opt[2] == 1
        dv = mat[:,1]
        ev = mat[:,2][1:end-1] # TODO: verify
        return Bidiagonal(dv, ev, :U)
    else
        error("Uncaught case")
    end
end
# function spdiags(dv, v, m, n)
#     Bidiagonal(dv, ev::V, uplo::Symbol)
# end
