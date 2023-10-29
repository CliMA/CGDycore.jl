  mutable struct IS 
    #PetscLayout  map;
    max::Int
    min::Int
    data
    total::Array{Int,1}
    nonlocal::Array{Int,1}
    local_offset::Int
    complement::IS          /* IS wrapping nonlocal indices. */
    IS() = (x = new(); x.complement = x)
    IS(max,
       min,
       data,
       toatal,
       nonlocal,
       local_offset,
       complement) = new(
    max,
    min,
    data,
    toatal,
    nonlocal,
    local_offset,
    complement)
  end


  function CreateStride(len::Int, pStart::Int, step::Int)
    is = IS()
    is.total=zeros(Int,len)
    for p = pStart + 1 : len
      is.total[p] = pStart + (p - 1) * step
    end
    return is
  end
