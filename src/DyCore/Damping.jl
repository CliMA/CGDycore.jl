function Damping!(F,W,Global)
  @inbounds for iz = Global.Grid.nz-1:-1:1
    zLoc=Global.Grid.z[iz+1]
    if zLoc>=Global.Grid.H-Global.Model.StrideDamp
      Damp = Global.Model.Relax*
        sin(0.5*pi*(1.0 - (Global.Grid.H - zLoc)/Global.Model.StrideDamp))^2
      F[iz]-=Damp*W[iz]
    else
      break
    end
  end
end

function Damping!(Fu,Fv,Fw,U,V,W,Ug,Vg,Global)
  @inbounds for iz = Global.Grid.nz-1:-1:1
    zLoc=Global.Grid.z[iz+1]
    if zLoc>=Global.Grid.H-Global.Model.StrideDamp
      Damp = Global.Model.Relax*
        sin(0.5*pi*(1.0 - (Global.Grid.H - zLoc)/Global.Model.StrideDamp))^2
      Fw[iz]-=Damp*W[iz]
      Fu[iz]-=Damp*(U[iz] - Ug[iz])
      Fv[iz]-=Damp*(V[iz] - Vg[iz])
    else
      break
    end
  end
end

function DampingKoeff!(K,CG,Global)
StrideDamp = Global.Model.StrideDamp
Relax = Global.Model.Relax
nz = Global.Grid.nz
z = Global.Grid.z
H = Global.Grid.H
K .= 0.0
for iz=nz-1:-1:1
  zLoc=z[iz+1]
  if zLoc>=H-StrideDamp
    K[iz] = -Relax *
      sin(0.5*pi*(1.0 - (H - zLoc)/StrideDamp))^2
  else
    break
  end
end
end
