function f = Div(Rhou,S,DS,dXdxIF,dXdxIC,JC)

M=size(Rhou,1);
N=size(Rhou,2);
OP=size(Rhou,3);
f=zeros(size(Rhou));

for i=1:N
  for j=1:M
    f(j,i,:) = reshape(DS*reshape(dXdxIC(j,i,:,1,1).*Rhou(j,i,:),...
      OP,1),1,1,OP)+...
      (dXdxIF(j,i+1,:,2,2).*S(j,i+1,:)-dXdxIF(j,i,:,2,2).*S(j,i,:));
    f(j,i,:)=f(j,i,:)./JC(j,i,:);
  end
end
end

