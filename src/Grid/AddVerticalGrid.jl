function AddVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
Grid.dz=H/nz;
Grid.zP[1]=Grid.dz/2;
Grid.H=H
for i=2:nz
  Grid.zP[i]=Grid.zP[i-1]+Grid.dz;
end
for i=2:nz+1
  Grid.z[i]=Grid.z[i-1]+Grid.dz;
end
end
