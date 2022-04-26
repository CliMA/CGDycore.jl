function vtkOutput(U,vtkGrid,CG,Global)
nz=Global.Grid.nz;
cOut=zeros(nz,CG.NumG,length(Global.Output.cNames));
for i=1:length(Global.Output.cNames)
  str = Global.Output.cNames[i]
  if str == "Rho"
      @views cOut[:,:,i]=U[:,:,Global.Model.RhoPos];
  elseif str == "u"
      @views cOut[:,:,i]=U[:,:,Global.Model.uPos];
  elseif str == "v"
      @views cOut[:,:,i]=U[:,:,Global.Model.vPos];
  elseif str == "w"
      BoundaryWOutput!(view(cOut,1,:,i),U,CG,Global)
      @views cOut[1,:,i] .= 0.5 .* (U[1,:,Global.Model.wPos] .+ cOut[1,:,i])
      @views cOut[2:nz-1,:,i] .= 0.5 .* (U[1:nz-2,:,Global.Model.wPos] .+ U[2:nz-1,:,Global.Model.wPos]);
      if nz>1
        @views cOut[nz,:,i] .=  0.5 .* U[nz-1,:,Global.Model.wPos];
      end
  elseif str == "Th"
      @views cOut[:,:,i].=U[:,:,Global.Model.ThPos]./U[:,:,Global.Model.RhoPos];
  elseif str == "RhoTh"
      @views cOut[:,:,i].=U[:,:,Global.Model.ThPos]
  elseif str == "TPrime"
      p=Pressure(U[:,:,Global.Model.ThPos],U[:,:,Global.Model.ThPos],U[:,:,Global.Model.ThPos],Global);
      @views cOut[:,:,i]=p./(Global.Rd*U[:,:,Global.Model.RhoPos])-Global.TBGrd;
  elseif str == "ThetaPrime"
      @views cOut[:,:,i]=U[:,:,Global.Model.ThPos]./U[:,:,Global.Model.RhoPos]-Global.ThetaBGrd;
  elseif str == "Pres"
      @views cOut[:,:,i]=Pressure(U[:,:,Global.Model.ThPos],U[:,:,Global.Model.ThPos],U[:,:,Global.Model.ThPos],Global);
  elseif str == "Vort"
      @views FVort2Vec!(cOut[:,:,i],U,CG,Global);
  end
end
vtk=vtkCG(cOut,CG,Global,vtkGrid,Global.Output.vtk);
return vtk
end
