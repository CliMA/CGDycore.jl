@muladd @inline function DerivativeX!(Dc::AbstractArray{FT,1},c::AbstractArray{FT,1},
  DX::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(DX,1)
  ny = size(DX,2)
  indDc = 0
  indj = 0
  indc = 0
  @inbounds for j = 1 : ny
    @inbounds for i = 1 : nx
      indDc += 1
      indc = indj
      @inbounds for l = 1 : nx
        indc += 1
        Dc[indDc] += DX[i,l] * c[indc]
      end
    end
    indj = indc
  end
end
@muladd @inline function DerivativeX!(Dc::AbstractArray{FT,2},c::AbstractArray{FT,2},
  DX::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(DX,1)
  ny = size(DX,2)
  indDc = 0
  indj = 0
  indc = 0
  @inbounds for j = 1 : ny
    @inbounds for i = 1 : nx
      indDc += 1
      indc = indj
      @inbounds for l = 1 : nx
        indc += 1
        @views @. Dc[indDc,:] += DX[i,l] * c[indc,:]
      end
    end
    indj = indc
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

@muladd @inline function DerivativeY!(Dc::AbstractArray{FT,1},c::AbstractArray{FT,1},
  DY::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(DY,1)
  ny = size(DY,2)
  indi = 0
  @inbounds for i = 1 : nx
    indi += 1
    indDc = indi
    @inbounds for j = 1 : ny
      indc = indi
      @inbounds for l = 1 : ny
        Dc[indDc] += DY[j,l] * c[indc]
        indc += nx
      end
      indDc += nx
    end
  end
end
  
@muladd @inline function DerivativeY!(Dc::AbstractArray{FT,2},c::AbstractArray{FT,2},
  DY::AbstractArray{FT,2}) where FT<:AbstractFloat
  nx = size(DY,1)
  ny = size(DY,2)
  indi = 0
  @inbounds for i = 1 : nx
    indi += 1
    indDc = indi
    @inbounds for j = 1 : ny
      indc = indi
      @inbounds for l = 1 : ny
        @views @. Dc[indDc,:] += DY[j,l] * c[indc,:]
        indc += nx
      end
      indDc += nx
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
