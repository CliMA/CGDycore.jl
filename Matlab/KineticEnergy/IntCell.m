function I= IntCell(c,JC,w)
%formula (34)
I=0;
M=size(c,1);
N=size(c,2);
OP=size(c,3);
for j=1:M
  for i=1:N
    I=I+w'*reshape(c(j,i,:).*JC(j,i,:),OP,1);
  end
end
end

