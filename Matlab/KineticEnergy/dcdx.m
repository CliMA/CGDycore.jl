function f = dcdx(c,DS,dXdxIC,JC)

M=size(c,1);
N=size(c,2);
OP=size(c,3);
f=zeros(size(c));

for i=1:N
  for j=1:M
    f(j,i,:) = dXdxIC(j,i,:,1,1).*reshape(DS*reshape(c(j,i,:),OP,1),1,1,OP)...
      ./JC(j,i,:);
  end
end

end
