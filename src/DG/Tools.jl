@muladd @inline function DerivativeX!(Dc::AbstractArray{FT,2},c::AbstractArray{FT,2},
  DX::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(c,1)
  ny = size(c,2)
  @inbounds for j = 1 : ny
    @inbounds for i = 1 : nx
      @inbounds for l = 1 : nx
        Dc[i,j] = Dc[i,j] + DX[i,l] * c[l,j]
      end
    end
  end 
end
@muladd @inline function DerivativeX!(Dc::AbstractArray{FT,3},c::AbstractArray{FT,3},
  DX::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(c,1)
  ny = size(c,2)
  @inbounds for j = 1 : ny
    @inbounds for i = 1 : nx
      @inbounds for l = 1 : nx
        @views @. Dc[i,j,:] = Dc[i,j,:] + DX[i,l] * c[l,j,:]
      end
    end
  end 
end
  
@muladd @inline function DerivativeY!(Dc::AbstractArray{FT,2},c::AbstractArray{FT,2},
  DY::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(c,1)
  ny = size(c,2)
  @inbounds for i = 1 : nx
    @inbounds for j = 1 : ny
      @inbounds for l = 1 : ny
        Dc[i,j] = Dc[i,j] + DY[j,l] * c[i,l]
      end
    end
  end
end

@muladd @inline function DerivativeY!(Dc::AbstractArray{FT,3},c::AbstractArray{FT,3},
  DY::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(c,1)
  ny = size(c,2)
  @inbounds for i = 1 : nx
    @inbounds for j = 1 : ny
      @inbounds for l = 1 : ny
        @views @. Dc[i,j,:] = Dc[i,j,:] + DY[j,l] * c[i,l,:]
      end
    end
  end
end
