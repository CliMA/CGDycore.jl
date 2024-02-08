function StiffMatrixFD(backend,FTB,FeF::HDivElement,FeT::ScalarElement,Grid,QuadOrd,Jacobi)
    QQ = FEMSei.QuadRule{FTB}(Grid.Type,backend,QuadOrd)
    Weights = QQ.Weights
    Points = QQ.Points
    fFRef  = zeros(FeT.Comp,FeF.DoF,length(Weights))
    fTRef  = zeros(FeT.Comp,FeT.DoF,length(Weights))

    for i = 1 : length(Weights)
        for iComp = 1 : FeT.Comp
            for iD = 1 : FeT.DoF
              fTRef[iComp,iD,i] = FeT.phi[iD,iComp](Points[i,1],Points[i,2])
            end
        end
        for iComp = 1 : FeT.Comp
            for iD = 1 : FeF.DoF
              fFRef[iComp,iD,i] = FeF.Divphi[iD,iComp](Points[i,1],Points[i,2])
            end
        end
    end
    RowInd = Int64[]
    ColInd = Int64[]
    Val = Float64[]

    for iF = 1 : Grid.NumFaces
      DivLoc = zeros(FeT.DoF,FeF.DoF)
        for i = 1 : length(Weights)
          DF, detJ = Jacobi(Grid.Type,QQ.Points[i,1],QQ.Points[i,2],Grid.Faces[iF], Grid)
          DivLoc = DivLoc + Weights[i]*(fTRef[:,:,i]'*fFRef[:,:,i])
        end
        for j = 1 : size(DivLoc,2)
          for i = 1 : size(DivLoc,1)
            push!(RowInd,FeT.Glob[i,iF])
            push!(ColInd,FeF.Glob[j,iF])
            push!(Val,DivLoc[i,j])
          end
        end
    end
    Div = sparse(RowInd, ColInd, Val)
    return Div
end
    
