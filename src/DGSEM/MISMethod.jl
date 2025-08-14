mutable struct MISStruct{FT<:AbstractFloat}
  nStage::Int
  beta::Array{FT, 2}
  alfa::Array{FT, 2}
  gamma::Array{FT, 2}
  d::Array{FT, 2}
end

function MISStruct{FT}() where FT<:AbstractFloat
  nStage = 0
  beta = zeros(FT,0,0)
  alfa = zeros(FT,0,0)
  gamma = zeros(FT,0,0)
  d = zeros(FT, 0)
  return MISStruct{FT}(
    nStage,
    beta,
    alfa,
    gamma,
    d,
  )
end

function MISStruct{FT}(Method) where FT<:AbstractFloat
  str = Method
  if str == "MISRK4"
    nStage = 4
	    beta = zeros(FT, nStage+1, nStage)
            alfa = zeros(FT, nStage+1, nStage)
            gamma = zeros(FT, nStage+1, nStage)
            d = zeros(FT, nStage+1, 1)
    	    beta[2, 1] = 0.38758444641450318
            beta[3, 1] = -2.5318448354142823E-002
            beta[3, 2] = 0.38668943087310403
            beta[4, 1] = 0.20899983523553325
            beta[4, 2] = -0.45856648476371231
            beta[4, 3] = 0.43423187573425748
            beta[5, 1] = -0.10048822195663100
            beta[5, 2] = -0.46186171956333327
            beta[5, 3] = 0.83045062122462809
            beta[5, 4] = 0.27014914900250392
    
            alfa[3, 2] = 0.52349249922385610
            alfa[4, 2] = 1.1683374366893629
            alfa[4, 3] = -0.75762080241712637
            alfa[5, 2] = -3.6477233846797109E-002
            alfa[5, 3] = 0.56936148730740477
            alfa[5, 4] = 0.47746263002599681
    
            gamma[3, 2] = 0.13145089796226542
            gamma[4, 2] = -0.36855857648747881
            gamma[4, 3] = 0.33159232636600550
            gamma[5, 2] = -6.5767130537473045E-002
            gamma[5, 3] = 4.0591093109036858E-002
            gamma[5, 4] = 6.4902111640806712E-002
    
		d[2] = beta[2, 1]
		d[3] = beta[3, 1] + beta[3, 2]
		d[4] = beta[4, 1] + beta[4, 2] + beta[4, 3]
		d[5] = beta[5, 1] + beta[5, 2] + beta[5, 3] + beta[5, 4]
  elseif str == "RK3"
    nStage = 3
    beta = zeros(FT, nStage+1, nStage)
    alfa = zeros(FT, nStage+1, nStage)
    gamma = zeros(FT, nStage+1, nStage)
    d = zeros(FT, nStage+1, 1)
    beta[2, 1] = 1/3
    beta[3, 2] = 1/2
    beta[4, 3] = 1
    d[2] = 1/3
    d[3] = 1/2
    d[4] = 1
#=            
   MIS.Method='RKJeb';
    MIS.nStage=3;
    MIS.A=zeros(MIS.nStage+1,MIS.nStage);
    MIS.G=zeros(MIS.nStage+1,MIS.nStage);
    MIS.D=zeros(MIS.nStage+1,MIS.nStage);
    MIS.d=zeros(MIS.nStage+1);
    MIS.A(2,1)=2.0492941060709863e-001;
    MIS.A(3,1)=-4.5477553356788974e-001;
    MIS.A(3,2)=9.5613538239378981e-001;
    MIS.A(4,1)=-3.5970281266252929e-002;
    MIS.A(4,2)=-1.5363649484946584e-001;
    MIS.A(4,3)=7.0259062712330234e-001;

    MIS.G(3,2)=-8.2176071248067006e-001;
    MIS.G(4,2)=-3.8080670922635063e-001;
    MIS.G(4,3)=4.5653105107801978e-001;

    MIS.D(3,2)=7.0302371060435331e-001;
    MIS.D(4,2)=4.2492220536139252e-001;
    MIS.D(4,3)=5.4545718243573982e-001;
=#   
   end
  return MISStruct{FT}(
    nStage,
    beta,
    alfa,
    gamma,
    d)
end

