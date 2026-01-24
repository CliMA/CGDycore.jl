mutable struct DiagonalMatrix{FT<:AbstractFloat} 
  m::Int
  n::Int
  D::Array{FT,1}    
end
mutable struct BlockDiagonalMatrix{FT<:AbstractFloat} 
  m::Int
  n::Int
  Index::Array{Array{Int,1},1}
  D::Array{SparseArrays.UMFPACK.UmfpackLU,1}    
end

function DiagonalMatrix(FT,FE)
  m = FE.M.m
  n = FE.M.n
  D = zeros(FT,m)
  for i = 1 : n
    D[i] = FE.M[i,i]
  end   
  return DiagonalMatrix{FT}(
    m,
    n,
    D
  )  
end  

function ldiv!(D::DiagonalMatrix,b)
  @. b /= D.D
end  


function BlockDiagonalDualMatrix(FT,FE,Grid)

  m = FE.M.m
  n = FE.M.n
  NumNodesN = 0 
  for iN = 1 : Grid.NumNodes
    if Grid.Nodes[iN].Type == 'N'
      NumNodesN += 1  
    end  
  end
  D = Array{SparseArrays.UMFPACK.UmfpackLU,1}(undef,NumNodesN)
  IndexGlob = Array{Array{Int, 1}, 1}(undef, NumNodesN)
  iM = 1
  for iN = 1 : Grid.NumNodes
    if Grid.Nodes[iN].Type == 'N'
      Index = []
      for iF in Grid.Nodes[iN].F
        for iDoF = 1 : FE.DoF
          push!(Index,FE.Glob[iDoF,iF])
        end
      end
      unique!(Index)
      IndexGlob[iM] = Index
      D[iM] = lu(FE.M[Index,Index])
      iM += 1
    end
  end
  return BlockDiagonalMatrix{FT}(
    m,
    n,
    IndexGlob,
    D
  )
end

function ldiv!(D::BlockDiagonalMatrix,b)
  @inbounds for iB = 1 : size(D.D,1)
    b[D.Index[iB]] = ldiv!(D.D[iB],b[D.Index[iB]])
  end   
end  
