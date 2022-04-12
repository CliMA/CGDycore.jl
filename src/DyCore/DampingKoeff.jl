function DampingKoeff!(K,CG,Global)
StrideDamp = Global.Model.StrideDamp
Relax = Global.Model.Relax
nz = Global.Grid.nz
z = Global.Grid.z
H = Global.Grid.H
K .= 0.0
for iz=nz-1:-1:1
  zLoc=z[iz+1];
  if zLoc>=H-StrideDamp
    K[iz] = -Relax *
      sin(0.5*pi*(1.0 - (H - zLoc)/StrideDamp))^2;
  else
    break
  end
end
end
