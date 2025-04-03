abstract type AdaptGrid end

Base.@kwdef struct GalChen <: AdaptGrid end

function (::GalChen)(H)
  @inline function AdaptHeight(zRef,zs)
    z = zRef + (H - zRef) * zs / H
    dzdzs = (H - zRef) / H
    dzdzRef = eltype(zRef)(1) - zs / H
    return z,dzdzref,dzdzs
  end
  return AdaptHeight
end

Base.@kwdef struct Sleve{T} <: AdaptGrid 
  etaH::T = .7
  s::T = 8/10

end

function (F::Sleve)(H)
  (;s,etaH) = F
  @inline function AdaptHeight(zRef,zs)
    eta = zRef / H
    if eta <= etaH
      z = eta * H + zs * sinh((etaH - eta) / s / etaH) / sinh(eltype(zRef)(1) / s) 
      dzdzRef = (H - zs * cosh((etaH - eta) / s / etaH) / sinh(eltype(zRef)(1) / s)) / H
      dzdzs = sinh((etaH - eta) / s / etaH) / sinh(eltype(zRef)(1) / s)
    else
      z = eta * H
      dzdzRef = eltype(zRef)(1)
      dzdzs = eltype(zRef)(0)
    end  
    return z,dzdzRef,dzdzs
  end
  return AdaptHeight
end

function AdaptGrid(FT,Type,H)

  if Type == "GalChen"
    AdaptGridFunction = Grids.GalChen()(H)
  elseif Type == "Sleve"
    AdaptGridFunction = Grids.Sleve{FT}()(H)
  end
  return AdaptGridFunction
end
