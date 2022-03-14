function vtkWriteHex(filename,X,Connectivity,c,cNames)

precision = '2';
c(abs(c)<=1.e-20)=0;
fid=fopen(filename, 'w');
fprintf(fid, '# vtk DataFile Version 2.0 \n');
fprintf(fid, 'Unstructured Grid Example \n');
fprintf(fid, 'ASCII \n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid, 'POINTS %d float \n',size(X,2));

spec = ['%0.', precision, 'f '];
fprintf(fid, ' %d %d %d \n', X);
%fprintf(fid, '\n');
fprintf(fid, 'CELLS %d %d \n',size(Connectivity,2),9*size(Connectivity,2));
fprintf(fid,'8 %d %d %d %d %d %d %d %d\n',Connectivity-1);
fprintf(fid, 'CELL_TYPES %d \n',size(Connectivity,2));
fprintf(fid,'%d  \n',12*ones(size(Connectivity,2),1));
fprintf(fid, 'Cell_DATA %d \n',size(Connectivity,2));
for i=1:size(c,2)
  fprintf(fid, 'SCALARS %s double 1 \n',cNames(i).s);
  fprintf(fid, 'LOOKUP_TABLE default \n');
  fprintf(fid,'%6e %6e %6e %6e \n',c(:,i));
end
fclose(fid);
end
