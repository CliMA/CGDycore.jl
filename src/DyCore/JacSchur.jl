mutable struct JStruct
    JRhoW::Array{Float64, 3}
    JWTh::Array{Float64, 3}
    JWRho::Array{Float64, 3}
    JThW::Array{Float64, 3}
    JTrW::Array{Float64, 4}
    JWW::Array{Float64, 3}
    tri::Array{Float64, 3}
    sw::Array{Float64, 2}
    CompTri::Bool
    CacheCol1::Array{Float64, 1}
    CacheCol2::Array{Float64, 1}
    CacheCol3::Array{Float64, 1}
end

function JStruct()
  JRhoW=zeros(0,0,0)
  JWTh=zeros(0,0,0)
  JWRho=zeros(0,0,0)
  JThW=zeros(0,0,0)
  JTrW=zeros(0,0,0,0)
  JWW=zeros(0,0,0)
  tri=zeros(0,0,0)
  sw=zeros(0,0)
  CompTri=false
  CacheCol1=zeros(0)
  CacheCol2=zeros(0)
  CacheCol3=zeros(0)
  return JStruct(
    JRhoW,
    JWTh,
    JWRho,
    JThW,
    JTrW,
    JWW,
    tri,
    sw,
    CompTri,
    CacheCol1,
    CacheCol2,
    CacheCol3,
  )
end  

function JStruct(NumG,nz,NumTr)
  JRhoW=zeros(2,nz,NumG)
  JWTh=zeros(2,nz,NumG)
  JWRho=zeros(2,nz,NumG)
  JThW=zeros(2,nz,NumG)
  JTrW=zeros(2,nz,NumG,NumTr)
  JWW=zeros(1,nz,NumG)
  tri=zeros(3,nz,NumG)
  sw=zeros(nz,NumG)
  CompTri=false
  CacheCol1=zeros(nz)
  CacheCol2=zeros(nz)
  CacheCol3=zeros(nz)
  return JStruct(
    JRhoW,
    JWTh,
    JWRho,
    JThW,
    JTrW,
    JWW,
    tri,
    sw,
    CompTri,
    CacheCol1,
    CacheCol2,
    CacheCol3,
  )
end

function JacSchur!(J,U,CG,Global,Param)
  (;  RhoPos,
      uPos,
      vPos,
      wPos,
      ThPos,
      NumV,
      NumTr) = Global.Model
  nz=Global.Grid.nz;
  NF=Global.Grid.NumFaces
  nCol=size(U,2);
  nJ=nCol*nz;

  D = J.CacheCol2
  Dp = J.CacheCol2
  Dm = J.CacheCol3
  dPdTh = J.CacheCol1
  K = J.CacheCol1

  @inbounds for iC=1:nCol
    @views Pres = Global.Cache.PresG[:,iC]
    @views Rho = U[:,iC,RhoPos]
    @views Th = U[:,iC,ThPos]
    @views Tr = U[:,iC,NumV+1:end]
    @views dz = Global.Metric.dz[:,iC]

    @views @. D[1:nz-1] = -(Rho[1:nz-1] + Rho[2:nz]) / (dz[1:nz-1] + dz[2:nz])
    @views @. J.JRhoW[1,:,iC] = D
    @views @. J.JRhoW[2,1:nz-1,iC] = -D[1:nz-1]

    dPresdTh!(dPdTh, Th, Rho, Tr, Global);
    @views @. Dp[1:nz-1] = dPdTh[2:nz] /  (0.25 * (Rho[1:nz-1] + Rho[2:nz])) / (dz[1:nz-1] + dz[2:nz])
    @views @. Dm[1:nz-1] = dPdTh[1:nz-1] / (0.25 * (Rho[1:nz-1] + Rho[2:nz])) / (dz[1:nz-1] + dz[2:nz])
    @views @. J.JWTh[1,2:nz,iC] = -Dp[1:nz-1]
    @views @. J.JWTh[2,:,iC] = Dm

    @views @. D[1:nz-1] = 0.5*(Pres[2:nz] - 2.0 * Pres[1:nz-1]) / (dz[1:nz-1] + dz[2:nz]) /
      (0.5*(Rho[1:nz-1] + Rho[2:nz]))^2;
    @views @. J.JWRho[1,2:nz,iC] = D[1:nz-1] 
    @views @. J.JWRho[2,1:nz,iC] = D 

    @views @. D[1:nz-1] = -0.5*(Th[1:nz-1] + Th[2:nz]) 
    if Global.Model.Thermo == "TotalEnergy" 
      @views @. D[1:nz-1] -= 0.5*(Pres[1:nz-1] + Pres[2:nz]) 
    end  
    @views @. J.JThW[1,:,iC] = D / dz
    @views @. J.JThW[2,1:nz-1,iC] = -D[1:nz-1] / dz[1:nz-1]

    @inbounds for iT = 1 : NumTr
      @views @. D[1:nz-1] = -0.5*(Tr[1:nz-1,iT] + Tr[2:nz,iT]) 
      @views @. J.JTrW[1,:,iC,iT] = D / dz
      @views @. J.JTrW[2,1:nz-1,iC,iT] = -D[1:nz-1] / dz[1:nz-1]
    end    

    if Global.Model.Damping
      @views DampingKoeff!(J.JWW[1,:,iC],CG,Global)
    end
  end
end

