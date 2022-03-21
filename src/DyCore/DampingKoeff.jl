function DampingKoeff(CG,Param)
K=zeros(CG.NumG,Param.Grid.nz);
H = Param.H
StrideDamp = Param.StrideDamp
Relax = Param.Relax
nz = Param.Grid.nz
z = Param.Grid.z
H = Param.H
for iz=nz-1:-1:1
  zLoc=z[iz+1];
  if zLoc>=H-StrideDamp
    K[:,iz] .= -Relax*
      sin(0.5*pi*(1.0 - (H - zLoc)/StrideDamp))^2;
  else
    break
  end
end
return K
end
