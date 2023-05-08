@muladd @inline function DerivativeX!(Dc::AbstractArray{Float64,2},c::AbstractArray{Float64,2},
  DX::AbstractArray{Float64,2})
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
@muladd @inline function DerivativeX!(Dc::AbstractArray{Float64,3},c::AbstractArray{Float64,3},
  DX::AbstractArray{Float64,2})
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
  
@muladd @inline function DerivativeY!(Dc::AbstractArray{Float64,2},c::AbstractArray{Float64,2},
  DY::AbstractArray{Float64,2})
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

@muladd @inline function DerivativeY!(Dc::AbstractArray{Float64,3},c::AbstractArray{Float64,3},
  DY::AbstractArray{Float64,2})
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
