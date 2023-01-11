function udwdx(u,wF,Rho,DS,dXdxIF,JC)
# formula (40)
  NF = size(wF,2)
  M = size(wF,1)
  OP = size(wF,3)
  f = zeros(size(wF))
  uAver = zeros(OP)
  tt = zeros(OP)
  for i = 2:NF-1
    for j=1:M
      @views @. uAver = (u[j,i-1,:] * Rho[j,i-1,:] * JC[j,i-1,:] +
        u[j,i,:] * Rho[j,i,:] * JC[j,i,:]) /
        (Rho[j,i-1,:] * JC[j,i-1,:] + Rho[j,i,:] * JC[j,i,:])
      @views tt = DS*reshape(wF[j,i,:],OP)
      @views @. f[j,i,:] = uAver * dXdxIF[j,i,:,1,1] * tt
      @views @. f[j,i,:] = 2.0 * f[j,i,:] / (JC[j,i-1,:]+JC[j,i,:])
    end
  end
  i=1
  for j=1:M
    @views @. uAver=(u[j,i,:] * Rho[j,i,:] * JC[j,i,:]) / (Rho[j,i,:] * JC[j,i,:])
    @views tt = DS*reshape(wF[j,i,:],OP)
    @views @. f[j,i,:] = uAver * dXdxIF[j,i,:,1,1] * tt
    @views @. f[j,i,:] = f[j,i,:] / (JC[j,i,:])
  end
  return f
end

