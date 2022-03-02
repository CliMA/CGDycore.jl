function Grid=HexagonGrid(nx,nz,lenx,lenz,Orient)
[Points,Faces,Cells]=Grid2DCFP(nx,nz,'hex',lenx,lenz);
Grid=CFPToGrid(Points,Faces,Cells,Orient);
end