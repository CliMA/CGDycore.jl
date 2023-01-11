function  wdwdx(wF,DS,dXdxIF,JC)
# formula (40)
  NF = size(wF,2)
  M = size(wF,1)
  OP = size(wF,3)
  N = NF-1
  f = zeros(M,N,OP)
  wdwLdx = zeros(OP)
  wdwRdx = zeros(OP)
  tt = zeros(OP)
  for i = 1:N
    for j = 1:M
      tt = DS*reshape(wF[j,i,:],OP)  
      @views @. wdwLdx = wF[j,i,:] * dXdxIF[j,i,:,1,1] * tt
      tt = DS*reshape(wF[j,i+1,:],OP)  
      @views @. wdwRdx = wF[j,i+1,:] * dXdxIF[j,i+1,:,1,1] * tt
      @views @. f[j,i,:] = 0.5 * (wdwLdx + wdwRdx) / JC[j,i,:]
    end
  end
  return f
end

