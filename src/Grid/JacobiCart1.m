function [D,J,X]=JacobiCart1(ksi,F,Grid)

x1=Grid.Nodes(F.N(1)).P(1);
y1=Grid.Nodes(F.N(1)).P(2);
x2=Grid.Nodes(F.N(2)).P(1);
y2=Grid.Nodes(F.N(2)).P(2);
x3=Grid.Nodes(F.N(3)).P(1);
y3=Grid.Nodes(F.N(3)).P(2);
x4=Grid.Nodes(F.N(4)).P(1);
y4=Grid.Nodes(F.N(4)).P(2);
xBar=ksi(1);
yBar=ksi(2);

J=[ (x1 - x2 + x3 - x4)*yBar + x2 - x1, (x1 - x2 + x3 - x4)*xBar + x4 - x1
     (y1 - y2 + y3 - y4)*yBar + y2 - y1, (y1 - y2 + y3 - y4)*xBar + y4 - y1];
J=[J
   0 0];
D=J(1,1)*J(2,2)-J(1,2)*J(2,1); 
if nargout > 2
  X=zeros(1,3);
  X(1)=Grid.Nodes(F.N(1)).P(1)...
    +(Grid.Nodes(F.N(2)).P(1)-Grid.Nodes(F.N(1)).P(1))*ksi(1)...
    +(Grid.Nodes(F.N(4)).P(1)-Grid.Nodes(F.N(1)).P(1))*ksi(2)...
    +(Grid.Nodes(F.N(3)).P(1)-Grid.Nodes(F.N(4)).P(1)-Grid.Nodes(F.N(2)).P(1)+Grid.Nodes(F.N(1)).P(1))*ksi(1)*ksi(2);
  X(2)=Grid.Nodes(F.N(1)).P(2)...
    +(Grid.Nodes(F.N(2)).P(2)-Grid.Nodes(F.N(1)).P(2))*ksi(1)...
    +(Grid.Nodes(F.N(4)).P(2)-Grid.Nodes(F.N(1)).P(2))*ksi(2)...
    +(Grid.Nodes(F.N(3)).P(2)-Grid.Nodes(F.N(4)).P(2)-Grid.Nodes(F.N(2)).P(2)+Grid.Nodes(F.N(1)).P(2))*ksi(1)*ksi(2);
  X(3)=Grid.Nodes(F.N(1)).P(3)...
    +(Grid.Nodes(F.N(2)).P(3)-Grid.Nodes(F.N(1)).P(3))*ksi(1)...
    +(Grid.Nodes(F.N(4)).P(3)-Grid.Nodes(F.N(1)).P(3))*ksi(2)...
    +(Grid.Nodes(F.N(3)).P(3)-Grid.Nodes(F.N(4)).P(3)-Grid.Nodes(F.N(2)).P(3)+Grid.Nodes(F.N(1)).P(3))*ksi(1)*ksi(2);
end
end
