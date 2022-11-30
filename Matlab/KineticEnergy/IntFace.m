function I= IntFace(cF,JC,w)
%formula (35)

I=0;
M=size(cF,1);
NF=size(cF,2);
OP=size(cF,3);
for j=1:M
  for i=2:NF-1
    I=I+w'*reshape(0.5*cF(j,i,:).*(JC(j,i-1,:)+JC(j,i,:)),OP,1);
  end
  I=I+w'*reshape(0.5*cF(j,1,:).*JC(j,1,:),OP,1);
  I=I+w'*reshape(0.5*cF(j,NF,:).*JC(j,NF-1,:),OP,1);
end
end

