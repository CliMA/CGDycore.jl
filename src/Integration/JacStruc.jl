struct JacStruct{A3,A4,A5,A6}
    JRhoW::A3
    JWRho::A3
    JThW::A3
    JWTh::A3
    JWW::A4
    tri::A5
    rRho::A6
    rTh::A6
    rw::A6
    sw::A6
end
function JacStruct(NumNodes,nz)

  JRhoW=zeros(2,NumNodes*nz)
  JWRho=zeros(2,NumNodes*nz)
  JThW=zeros(2,NumNodes*nz)
  JWTh=zeros(2,NumNodes*nz)
  JWW=zeros(1,NumNodes*nz)
  tri=zeros(3,NumNodes*nz)
  rRho=zeros(NumNodes*nz,1)
  rTh=zeros(NumNodes*nz,1)
  rw=zeros(NumNodes*nz,1)
  sw=zeros(NumNodes*nz,1)
  A3 = typeof(JRhoW)
  A4 = typeof(JWW)
  A5 = typeof(tri)
  A6 = typeof(sw)
  return JacStruct{A3,A4,A5,A6}(
        JRhoW,
        JWRho,
        JThW,
        JWTh,
        JWW,
        tri,
        rRho,
        rTh,
        rw,
        sw,
    )
end  

