function [cE] = DSS(c,JC)
N=size(c,2);
M=size(c,1);
OP=size(c,3);
cE=c;
for i=2:M
  cE(i-1,:,OP)=(c(i-1,:,OP).*JC(i-1,:,OP)+c(i,:,1).*JC(i,:,1))...
    ./(JC(i-1,:,OP)+JC(i,:,1));
  cE(i,:,1)=cE(i-1,:,OP);
end
cE(M,:,OP)=(c(M,:,OP).*JC(M,:,OP)+c(1,:,1).*JC(1,:,1))...
    ./(JC(M,:,OP)+JC(1,:,1));
cE(1,:,1)=cE(M,:,OP);
end

