using Printf
function vtkWriteHex(filename,X,Connectivity,c,cNames)

precision = "2";
c .= c .* (abs.(c) .> 1.e-20);
open(filename, "w") do fid
  print(fid, "# vtk DataFile Version 2.0 \n");
  print(fid, "Unstructured Grid Example \n");
  print(fid, "ASCII \n");
  print(fid, "DATASET UNSTRUCTURED_GRID \n");
  @printf(fid, "POINTS %d float \n",size(X,2));

  for i in 1:size(X,2)
    @printf(fid, " %6e %6e %6e \n", X[:,i]...);
  end
  @printf(fid, "\n");
  @printf(fid, "CELLS %d %d \n",size(Connectivity,2),9*size(Connectivity,2));
  for i in 1:size(Connectivity, 2)
    @printf(fid, "8 %d %d %d %d %d %d %d %d\n",(Connectivity[:,i] .- 1)...);
  end
  @printf(fid, "CELL_TYPES %d \n",size(Connectivity,2));
  for i in 1:size(Connectivity,2)
    @printf(fid, "%d  \n",12);
  end
  @printf(fid, "Cell_DATA %d \n",size(Connectivity,2));
  for i=1:size(c,2)
    @printf(fid, "SCALARS %s double 1 \n",cNames[i]);
    print(fid, "LOOKUP_TABLE default \n");
    for j=1:4:size(c,1)-3
      @printf(fid, "%6e %6e %6e %6e \n",c[j:j+3,i]...);
    end
  end
end
end

function vtkWriteQuad(filename,X,Connectivity,c,cNames)

precision = "2";
c .= c .* (abs.(c) .> 1.e-20);
open(filename, "w") do fid
  print(fid, "# vtk DataFile Version 2.0 \n");
  print(fid, "Unstructured Grid Example \n");
  print(fid, "ASCII \n");
  print(fid, "DATASET UNSTRUCTURED_GRID \n");
  @printf(fid, "POINTS %d float \n",size(X,2));

  for i in 1:size(X,2)
    @printf(fid, " %6e %6e %6e \n", X[:,i]...);
  end
  @printf(fid, "\n");
  @printf(fid, "CELLS %d %d \n",size(Connectivity,2),5*size(Connectivity,2));
  for i in 1:size(Connectivity, 2)
    @printf(fid, "4 %d %d %d %d\n",(Connectivity[:,i] .- 1)...);
  end
  @printf(fid, "CELL_TYPES %d \n",size(Connectivity,2));
  for i in 1:size(Connectivity,2)
    @printf(fid, "%d  \n",9);
  end
  @printf(fid, "Cell_DATA %d \n",size(Connectivity,2));
  for i=1:size(c,2)
    @printf(fid, "SCALARS %s double 1 \n",cNames[i]);
    print(fid, "LOOKUP_TABLE default \n");
    for j=1:4:size(c,1)-3
      @printf(fid, "%6e %6e %6e %6e \n",c[j:j+3,i]...);
    end
  end
end

end
