mutable struct ColumnElementStruct{N,T}
  Rho :: MArray{N,T}
  u::MArray{N,T}
  v::MArray{N,T}
  w::MArray{N,T}
  Theta::MArray{N,T}
  Tr::MArray{N,T}
end

function ColumnElement(NPOLY,NZ,T)
  Rho=@MArray zeros(T,NPOLY,NZ)
  u=@MArray zeros(T,NPOLY,NZ)
  v=@MArray zeros(T,NPOLY,NZ)
  w=@MArray zeros(T,NPOLY,NZ)
  Theta=@MArray zeros(T,NPOLY,NZ)
  Tr=@MArray zeros(T,NPOLY,NZ)

  return ColumnElementStruct(
    Rho,
    v,
    u,
    w,
    Theta,
    Tr,
  )
end  

function test1!(V1::ColumnElementStruct,V2::ColumnElementStruct,NPOLY,NZ)
  
  @SArray a = zeros(Float64,NPOLY,NZ)
  @show typeof(a)
  @. a = V2.Rho * V2.Rho
  @. V1.Rho += a
end  
