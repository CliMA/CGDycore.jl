function []=PlotFaceKiteGrid(Grid,fig)
figure(fig)
clf(fig)
hold on;
for iF=1:size(Grid.Faces,2)
  for iK=1:size(Grid.Faces(iF).Kite.Faces,2)
    X=zeros(size(Grid.Faces(iF).Kite.Faces(iK).N,2),3);
    for i=1:size(Grid.Faces(iF).Kite.Faces(iK).N,2)
      X(i,:)=Grid.Faces(iF).Kite.Nodes(Grid.Faces(iF).Kite.Faces(iK).N(i)).P;
    end
    fill3(X(:,1),X(:,2),X(:,3),0);
  end
end
hold off
end
