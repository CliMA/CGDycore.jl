function DampingKoeff(CG,Param)
K=zeros(CG.NumG,Param.Grid.nz);
for iz=Param.Grid.nz-1:-1:1
  zLoc=Param.Grid.z[iz+1];
  if zLoc>=Param.H-Param.StrideDamp
    K[:,iz] .= -Param.Relax*
      sin(0.5*pi*(1.0 - (Param.H - zLoc)/Param.StrideDamp))^2;
  else
    break
  end
end
return K
end
