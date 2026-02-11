abstract type DampingValue end

Base.@kwdef struct DampingW <: DampingValue end

function (::DampingW)(H,StrideDamp,Relax,uPos,vPos,wPos,::Examples.VelocityS,::Grids.SphericalGrid)
  @inline function Damping(X,U)
    FT = eltype(X)
    Rad = sqrt(X[1]^2 + X[2]^2 + X[3]^2)
    z = Rad - FT(P.RadEarth)
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

function (::DampingW)(H,StrideDamp,Relax,uPos,vPos,wPos,::Examples.VelocityC,::Grids.SphericalGrid)
  @inline function Damping(X,U)
    FT = eltype(X)
    Rad2 = X[1]^2 + X[2]^2 + X[3]^2
    Rad = sqrt(Rad2)
    z = Rad - P.RadEarth
    if z>=H-StrideDamp
      Damp = Relax *
        sin(FT(0.5) * pi * (FT(1) - (H - z)/StrideDamp))^2
      Un = (X[1] * U[uPos] + X[2] * U[vPos] + X[2] * U[wPos]) / Rad2  
      Fu = -Damp * Un * X[1]
      Fv = -Damp * Un * X[2]
      Fw = -Damp * Un * X[3]
    else
      Fu = FT(0)
      Fv = FT(0)
      Fw = FT(0)
    end
    return Fu,Fv,Fw
  end
  return Damping
end

function (::DampingW)(H,StrideDamp,Relax,uPos,vPos,wPos,::Examples.VelocityC,::Grids.CartesianGrid)
  @inline function Damping(X,U)
    FT = eltype(X)
    z = X[3] 
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


