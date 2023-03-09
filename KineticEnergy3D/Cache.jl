mutable struct Cache
  pC::Array{Float64, 6}
  pF::Array{Float64, 6}
  UUF::Array{Float64, 7}
  FUUF::Array{Float64, 7}
  KinF::Array{Float64, 5}
  KinC::Array{Float64, 5}
  CF::Array{Float64, 7}
  Column::Array{Float64, 5}
  Block::Array{Float64, 4}
  BlockXY::Array{Float64, 3}
end

function Cache(Nx,Ny,Nz,OrdPolyX,OrdPolyY,OrdPolyZ)
  pC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ)
  pF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1)
  UUF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,4)
  FUUF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,4)
  KinF = zeros(Nx,Ny,Nz+1,OrdPolyX+1,OrdPolyY+1)
  KinC = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1)
  CF = zeros(Nx,Ny,Nz,OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,5)
  Column = zeros(OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,Nz,3)
  Block = zeros(OrdPolyX+1,OrdPolyY+1,OrdPolyZ+1,6)
  BlockXY = zeros(OrdPolyX+1,OrdPolyY+1,3)

  return Cache(
    pC,
    pF,
    UUF,
    FUUF,
    KinF,
    KinC,
    CF,
    Column,
    Block,
    BlockXY,
  )
end  
