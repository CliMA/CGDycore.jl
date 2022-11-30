function f = wdwdx(wF,DS,dXdxIF,JC)
%formula (40)
NF=size(wF,2);
M=size(wF,1);
OP=size(wF,3);
N=NF-1;
f=zeros(M,N,OP);
for i=1:N
  for j=1:M
    wdwLdx=wF(j,i,:).*dXdxIF(j,i,:,1,1).*reshape(DS*reshape(wF(j,i,:),OP,1),1,1,OP);
    wdwRdx=wF(j,i+1,:).*dXdxIF(j,i+1,:,1,1).*reshape(DS*reshape(wF(j,i+1,:),OP,1),1,1,OP);
    f(j,i,:) = 0.5*(wdwLdx+wdwRdx)./JC(j,i,:);
  end
end
end

