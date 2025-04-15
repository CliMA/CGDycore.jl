using DynamicPolynomials
using FastGaussQuadrature
using LinearAlgebra
function Polynomial_k(k,x)
  DoF::Int = (k + 1) *(k + 2) / 2
  phi = Array{Polynomial,1}(undef,DoF)
  iDoF = 1
  for i = 0 : k
    for j = 0 : i
      phi[iDoF] = x[1][1]^(i-j) * x[1][2]^j + 0.0
      iDoF += 1
    end
  end
  return phi
end

function QuadPoints(Method)

  N1 = [0.0 -1.0]
  N2 = [1.0 1.0]
  N3 = [-1.0 0.0]
  if Method  == "Kubatko1"
    k = 1  
    n = 6
    ksi = zeros(2,n)
    N = zeros(2,n)
    ksi[:,1] = [-0.577350269189626  -1]
    ksi[:,2] = [0.577350269189626 -1]
    ksi[:,3] = [-1 -0.577350269189626]
    ksi[:,4] = [-1 0.577350269189626]
    ksi[:,5] = [-0.577350269189626 0.577350269189626]
    ksi[:,6] = [ 0.577350269189626 -0.577350269189626]

    N[:,1] = N1
    N[:,2] = N1
    N[:,3] = N3
    N[:,4] = N3
    N[:,5] = N2
    N[:,6] = N2

    w = ones(n) * 1/3
    wF = ones(n)
  elseif Method  == "Kubatko2"
    k = 2
    n = 10
    ksi = zeros(2,n)
    N = zeros(2,n)
    w = zeros(n)
    wF = zeros(n)

    ksi[:,1] = [-0.774596669241483 -1] 
    N[:,1] = N1
    ksi[:,2] = [0 -1] 
    N[:,2] = N1
    ksi[:,3] = [0.774596669241483 -1] 
    N[:,3] = N1
    ksi[:,4] = [-1 -0.774596669241483] 
    N[:,4] = N3
    ksi[:,5] = [-1 0] 
    N[:,5] = N3
    ksi[:,6] = [-1 0.774596669241483]
    N[:,6] = N3
    ksi[:,7] = [-0.774596669241483 0.774596669241483] 
    N[:,7] = N2
    ksi[:,8] = [0 0] 
    N[:,8] = N2
    ksi[:,9] = [0.774596669241483 -0.774596669241483] 
    N[:,9] = N2
    ksi[:,10] = [-0.333333333333333 -0.333333333333333]

    w[1] = 0.083333333333333 
    w[2] = 0.2 
    w[3] = 0.083333333333333 
    w[4] = 0.083333333333333 
    w[5] = 0.2
    w[6] = 0.083333333333333 
    w[7] = 0.083333333333333 
    w[8] = 0.2 
    w[9] = 0.083333333333333 
    w[10] = 0.9


    _, wGL = gausslegendre(k+1)
    wF[1:3] = wGL
    wF[4:6] = wGL
    wF[7:9] = wGL
    
  elseif Method == "Hicken1"
    k = 1
    n = 7
    ksi = zeros(2,n)
    N = zeros(2,n)
    w = zeros(n)
    wF = zeros(n)
    wFx1 = zeros(n)
    wFx1 = zeros(n)
    ksi[:,1] = [-1.0 -1.0]
    N[:,1] = [0.0 -1.0] + [-1.0 0.0]
    ksi[:,2] = [1.0 -1.0]
    N[:,2] = [0.0 -1.0] + [1.0 1.0]
    ksi[:,3] = [-1.0 1.0]
    N[:,3] = [1.0 1.0] + [-1.0 0.0]
    ksi[:,4] = [0.0 -1.0]
    N[:,4] = [0.0 -1.0]
    ksi[:,5] = [0.0 0.0]
    N[:,5] = [1.0 1.0]
    ksi[:,6] = [-1.0 0.0]
    N[:,6] = [-1.0 0.0]
    ksi[:,7] = [-1/3 -1/3]

    w[1] = 0.09999999999999999
    w[2] =0.09999999999999999
    w[3] =0.09999999999999999
    w[4] =0.26666666666666666
    w[5] =0.26666666666666666
    w[6] =0.26666666666666666
    w[7] =0.9000000000000002
    wF = zeros(n)
    wF[1] = 0.333333333333333 
    wF[2] = 0.333333333333333 
    wF[3] = 0.333333333333333 
    wF[4] = 1.333333333333333
    wF[5] = 1.333333333333333
    wF[6] = 1.333333333333333
  end
  wFx1 = wF .* N[1,:]
  wFx2 = wF .* N[2,:]
  s = @polyvar x[1:2]
  phi = Polynomial_k(k,s)
  nSt = size(phi,1)
  phiDx1 = Array{Polynomial,2}(undef,nSt,1)
  phiDx2 = Array{Polynomial,2}(undef,nSt,1)
  for i = 1 : nSt
    phiDx1[i] = differentiate(phi[i],x[1])  
    phiDx2[i] = differentiate(phi[i],x[2])  
  end  
  V = zeros(n,nSt)
  VDx1 = zeros(n,nSt)
  VDx2 = zeros(n,nSt)
  for j = 1 : nSt
    for i = 1 : n
      V[i,j] = phi[j](ksi[1,i],ksi[2,i])  
      VDx1[i,j] = phiDx1[j](ksi[1,i],ksi[2,i])  
      VDx2[i,j] = phiDx2[j](ksi[1,i],ksi[2,i])  
    end
  end  
  Rx1 = V' * diagm(w) * VDx1 + VDx1' * diagm(w) * V
  RFx1 = V' * diagm(wFx1) * V
  @show sum(abs.(RFx1-Rx1))
  Rx2 = V' * diagm(w) * VDx2 + VDx2' * diagm(w) * V
  RFx2 = V' * diagm(wFx2) * V
  @show sum(abs.(RFx2-Rx2))


  Q,R = qr(diagm(sqrt.(w)) * V)
  VO = V / R 
  VODx1 = VDx1 / R 
  VODx2 = VDx2 / R 

  Dx1 = 0.5 * (diagm(1.0 ./ w) + VO * VO') * diagm(wFx1) * (I - VO * VO' * diagm(w)) +
    VODx1 * V' * diagm(w)
  Dx2 = 0.5 * (diagm(1.0 ./ w) + VO * VO') * diagm(wFx2) * (I - VO * VO' * diagm(w)) +
    VODx2 * V' * diagm(w)

  return Dx1, Dx2

end

k = 1
Dx1K1, Dx2K1 = QuadPoints("Kubatko1")
Dx1K2, Dx2K2 = QuadPoints("Kubatko2")
Dx1H1, Dx2H1 = QuadPoints("Hicken1")
a=2
