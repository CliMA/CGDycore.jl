abstract type DampingValue end

Base.@kwdef struct DampingW <: DampingValue end

function (::DampingW)(H,StrideDamp,Relax,wPos)
  @inline function Damping(z,U)
    FT = eltype(z)
    if z>=H-StrideDamp
      Damp = Relax *
        sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
      Fw = -Damp * U[wPos]
    else
      Fw = FT(0)
    end
    return FT(0),FT(0),Fw
  end
  return Damping
end

