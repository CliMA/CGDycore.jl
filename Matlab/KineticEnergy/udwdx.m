function f = udwdx(u,wF,Rho,DS,dXdxIF,JC)
%formula (40)
NF=size(wF,2);
M=size(wF,1);
OP=size(wF,3);
f=zeros(size(wF));
for i = 2:NF-1
  for j=1:M
    uAver=(u(j,i-1,:).*Rho(j,i-1,:).*JC(j,i-1,:)...
      +u(j,i,:).*Rho(j,i,:).*JC(j,i,:))...
      ./(Rho(j,i-1,:).*JC(j,i-1,:)+Rho(j,i,:).*JC(j,i,:));
     f(j,i,:) = uAver...
       .*dXdxIF(j,i,:,1,1).*reshape(DS*reshape(wF(j,i,:),OP,1),1,1,OP);
     f(j,i,:) = 2.0*f(j,i,:)./(JC(j,i-1,:)+JC(j,i,:));
  end
end
i=1;
for j=1:M
  uAver=(...
    +u(j,i,:).*Rho(j,i,:).*JC(j,i,:))...
    ./(Rho(j,i,:).*JC(j,i,:));
  f(j,i,:) = uAver...
    .*dXdxIF(j,i,:,1,1).*reshape(DS*reshape(wF(j,i,:),OP,1),1,1,OP);
  f(j,i,:) = 2.0*f(j,i,:)./(JC(j,i,:));
end
end

