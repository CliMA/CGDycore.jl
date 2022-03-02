function [D,J,X]=JacobiTri(ksi,F,Grid)
dX1dksi1=(Grid.Nodes(F.N(2)).P(1)-Grid.Nodes(F.N(1)).P(1));
dX1dksi2=(Grid.Nodes(F.N(3)).P(1)-Grid.Nodes(F.N(1)).P(1));
dX2dksi1=(Grid.Nodes(F.N(2)).P(2)-Grid.Nodes(F.N(1)).P(2));
dX2dksi2=(Grid.Nodes(F.N(3)).P(2)-Grid.Nodes(F.N(1)).P(2));
dX3dksi1=(Grid.Nodes(F.N(2)).P(3)-Grid.Nodes(F.N(1)).P(3));
dX3dksi2=(Grid.Nodes(F.N(3)).P(3)-Grid.Nodes(F.N(1)).P(3));
J=zeros(3,2);

J(1,1)=dX1dksi1;
J(2,1)=dX2dksi1;
J(3,1)=dX3dksi1;
J(1,2)=dX1dksi2;
J(2,2)=dX2dksi2;
J(3,2)=dX3dksi2;

D=norm(cross(J(:,1),J(:,2)),2);
D=J(1,1)*J(2,2)-J(2,1)*J(1,2);
if nargout > 2
  X=zeros(1,3);
  X(1)=Grid.Nodes(F.N(1)).P(1)...
    +(Grid.Nodes(F.N(2)).P(1)-Grid.Nodes(F.N(1)).P(1))*ksi(1)...
    +(Grid.Nodes(F.N(3)).P(1)-Grid.Nodes(F.N(1)).P(1))*ksi(2);
    
  X(2)=Grid.Nodes(F.N(1)).P(2)...
    +(Grid.Nodes(F.N(2)).P(2)-Grid.Nodes(F.N(1)).P(2))*ksi(1)...
    +(Grid.Nodes(F.N(3)).P(2)-Grid.Nodes(F.N(1)).P(2))*ksi(2);
  X(3)=Grid.Nodes(F.N(1)).P(3)...
    +(Grid.Nodes(F.N(2)).P(3)-Grid.Nodes(F.N(1)).P(3))*ksi(1)...
    +(Grid.Nodes(F.N(3)).P(3)-Grid.Nodes(F.N(1)).P(3))*ksi(2);
end
end
