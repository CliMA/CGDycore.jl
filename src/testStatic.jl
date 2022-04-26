using StaticArrays
const NPOLY=5
const NZ=20

struct Cell
  DX11::MArray
end  

function Cell()
  DX11A=zeros(NPOLY,NPOLY)
  return(
  SMatrix(DX11A),
  )
end



