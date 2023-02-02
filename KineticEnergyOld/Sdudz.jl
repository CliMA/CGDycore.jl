function Sdudz(u,Rho,S,dXdxIF,JC)
#formula (40)
  NF = size(S,2)
  N = NF - 1
  M = size(S,1)
  f = zeros(size(u))
  uE = zeros(size(u,1),N+2,size(u,3))
  @views @. uE[:,2:N+1,:] = u[:,1:N,:]
  @views @. uE[:,1,:] = u[:,1,:]
  @views @. uE[:,N+2,:] = u[:,N,:]
  for i = 1 : N
    for j = 1 : M
      @views @. f[j,i,:] = 1.0 / Rho[j,i,:] *
      ((S[j,i+1,:] * dXdxIF[j,i+1,:,2,2] * (uE[j,i+2,:] + uE[j,i+1,:]) -
       S[j,i,:] * dXdxIF[j,i,:,2,2] * (uE[j,i+1,:] + uE[j,i,:])) /
      (2.0*JC[j,i,:]) -
      uE[j,i+1,:] * (S[j,i+1,:] * dXdxIF[j,i+1,:,2,2] - S[j,i,:] * dXdxIF[j,i,:,2,2]) /
      JC[j,i,:])
    end
  end
  return f
end

