function Grid=DGNodes(Grid,Trans,DG)
for iE=1:Grid.NumEdges
  Grid.Edges(iE).DGPoints=Trans(DG.xw,P1,P2,Grid);
end
end