function [cE] = DSSF(c,JC)
NF=size(c,2);
N=NF-1;
M=size(c,1);
OP=size(c,3);
cE=c;
JF=zeros(M,NF,OP);
JF(:,2:NF-1,:)=0.5*(JC(:,1:N-1,:)+JC(:,2:N,:));
JF(:,1,:)=JC(:,1,:);
JF(:,NF,:)=JC(:,N,:);
for j=2:M
  for i=1:NF
    cE(j-1,i,OP)=(c(j-1,i,OP)*JF(j-1,i,OP)+c(j,i,1)*JF(j,i,1))...
      /(JF(j-1,i,OP)+JF(j,i,1));
    cE(j,i,1)=cE(j-1,i,OP);
  end
end
for i=1:NF
  cE(M,i,OP)=(c(M,i,OP)*JF(M,i,OP)+c(1,i,1)*JF(1,i,1))...
    /(JF(M,i,OP)+JF(1,i,1));
  cE(1,i,1)=cE(M,i,OP);
end
end

