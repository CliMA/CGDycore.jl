function DSS!(cCG,CG,Grid)
  OPM1=CG.OrdPoly-1
  for iE = 1:Grid.NumEdges
    iF1=Grid.Edges[iE].F[1]  
    iF2=Grid.Edges[iE].F[2]  
    for i = 2:OPM1
      temp=cCG[1,i,:,iF1]  
      cCG[1,i,:,iF1]+=cCG[1,i,:,iF2]
      cCG[1,i,:,iF2]+=temp
    end
  end
  for iN = 1:Grid.NumNodes
    iF1=Grid.Nodes[iN].F[1]
    temp = 0.0
    for iF in Grid.Nodes[iN].F
      i=Grid.Nodes[iN].Ind1[iF]  
      j=Grid.Nodes[iN].Ind2[iF]  
      temp += cCG[i,j,:,iF]  
    end
    for iF in Grid.Nodes[iN].F
      i=Grid.Nodes[iN].Ind1[iF]  
      j=Grid.Nodes[iN].Ind2[iF]  
      cCG[i,j,:,iF] = temp  
    end
  end  

