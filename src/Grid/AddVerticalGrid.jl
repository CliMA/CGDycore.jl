function AddVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
@. Grid.dzeta = H/nz;
Grid.H=H
for i=2:nz+1
  Grid.z[i] = Grid.z[i-1] + Grid.dzeta[i-1]
end
for i=1:nz
  Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
end
end
