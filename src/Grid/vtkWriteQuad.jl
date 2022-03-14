function vtkWriteQuad(filename,X,Connectivity,c)

precision = '2';

fid=fopen(filename, 'w');
fprintf(fid, '# vtk DataFile Version 2.0 \n');
fprintf(fid, 'Unstructured Grid Example \n');
fprintf(fid, 'ASCII \n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid, 'POINTS %d float \n',size(X,2));

spec = ['%0.', precision, 'f '];
fprintf(fid, spec, X);
fprintf(fid, '\n');
fprintf(fid, 'CELLS %d %d \n',size(Connectivity,2),5*size(Connectivity,2));
fprintf(fid,'4 %d %d %d %d \n',Connectivity-1);
fprintf(fid, 'CELL_TYPES %d \n',size(Connectivity,2));
fprintf(fid,'%d  \n',9*ones(size(Connectivity,2),1));
fprintf(fid, 'Cell_DATA %d \n',size(Connectivity,2));
fprintf(fid, 'SCALARS scalars float 1 \n');
fprintf(fid, 'LOOKUP_TABLE default \n');
fprintf(fid,'%d %d %d %d \n',c);
fclose(fid);
end
