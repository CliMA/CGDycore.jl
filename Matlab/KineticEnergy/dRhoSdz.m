function f = dRhoSdz(S,dXdxIF,JC)

N=size(S,2)-1;
M=size(S,1);
f=zeros(size(S,1),N,size(S,3));

for i = 1:N
  for j=1:M
    f(j,i,:)=(dXdxIF(j,i+1,:,2,2).*S(j,i+1,:)-dXdxIF(j,i,:,2,2).*S(j,i,:))...
      ./JC(j,i,:);
  end
end
end

