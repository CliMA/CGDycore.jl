mutable struct MISStruct
  nStage::Int
  A::Array{Float64, 2}
  G::Array{Float64, 2}
  D::Array{Float64, 2}
  d::Array{Float64, 1}
end
function MISMethod()
  nStage=0
  A=zeros(0,0)
  G=zeros(0,0)
  D=zeros(0,0)
  d=zeros(0)
  return MISStruct(
    nStage,
    A,
    G,
    D,
    d,
  )
end

function MISMethod(Method)
  str = Method
  if str == "MISRK3"
    nStage = 3
    A = [0 0 0 0
         1.0 0 0 0
         0.0 1.0 0.0 0.0
         0.0 0.0 1.0 0.0]
    G=zeros(nStage+1,nStage)
    D=zeros(nStage+1,nStage)
    d = [0, 1.0/3.0, 1.0/2.0, 1.0]
  elseif str == "MISRKJeb"
    nStage=3
    A=zeros(nStage+1,nStage)
    G=zeros(nStage+1,nStage)
    D=zeros(nStage+1,nStage)
    d=zeros(nStage+1)
    A[2,1]=2.0492941060709863e-001
    A[3,1]=-4.5477553356788974e-001
    A[3,2]=9.5613538239378981e-001
    A[4,1]=-3.5970281266252929e-002
    A[4,2]=-1.5363649484946584e-001
    A[4,3]=7.0259062712330234e-001
    
    G[3,2]=-8.2176071248067006e-001
    G[4,2]=-3.8080670922635063e-001
    G[4,3]=4.5653105107801978e-001
    
    D[3,2]=7.0302371060435331e-001
    D[4,2]=4.2492220536139252e-001
    D[4,3]=5.4545718243573982e-001
    
    for i=2:nStage+1
      for j=1:i-1
        d[i]=d[i]+A[i,j]
      end
      for j=1:i-1
        A[i,j]=A[i,j]/d[i]
        G[i,j]=G[i,j]/d[i]
      end
    end
  end  
  return MISStruct(
    nStage,
    A,
    G,
    D,
    d,
  )
end
