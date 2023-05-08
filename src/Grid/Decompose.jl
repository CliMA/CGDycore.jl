function Decompose(Grid,OwnProc,NumProc)

NumFaces=Grid.NumFaces
LocalNumfaces=zeros(Int,NumProc)
LocalNumfaces.=NumFaces / NumProc
Rest = mod(NumFaces, NumProc)
for iP=1:Rest
  LocalNumfaces[iP]+=1
end  
CellLocalToGlobal = sum(LocalNumfaces[1:OwnProc-1])+1:sum(LocalNumfaces[1:OwnProc])
stop
    
end
