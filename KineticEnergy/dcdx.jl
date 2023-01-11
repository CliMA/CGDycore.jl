function dcdx(c,DS,dXdxIC,JC)

  M=size(c,1)
  N=size(c,2)
  OP=size(c,3)
  f=zeros(size(c))
  tt=zeros(OP)

  for i=1:N
    for j=1:M
      tt = DS*reshape(c[j,i,:],OP)  
      @views @. f[j,i,:] = dXdxIC[j,i,:,1,1] * tt / JC[j,i,:]
    end
  end
  return f
end

