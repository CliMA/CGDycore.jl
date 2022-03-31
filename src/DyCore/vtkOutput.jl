function vtkOutput(U,vtkGrid,CG,Param)
nz=Param.Grid.nz;
cOut=zeros(CG.NumG,nz,length(Param.cNames));
for i=1:length(Param.cNames)
  str = Param.cNames[i]
  if str == "Rho"
      @views cOut[:,:,i]=U[:,:,Param.RhoPos];
  elseif str == "u"
      @views cOut[:,:,i]=U[:,:,Param.uPos];
  elseif str == "v"
      @views cOut[:,:,i]=U[:,:,Param.vPos];
  elseif str == "w"
      BoundaryWOutput!(cOut[:,1,i],U,CG,Param)
      @views cOut[:,1,i] .= 0.5 .* (U[:,1,Param.wPos] .+ cOut[:,1,i])
      @views cOut[:,2:nz-1,i] .= 0.5 .* (U[:,1:nz-2,Param.wPos] .+ U[:,2:nz-1,Param.wPos]);
      @views cOut[:,nz,i] .=  0.5 .* U[:,nz-1,Param.wPos];
  elseif str == "Th"
      @views cOut[:,:,i].=U[:,:,Param.ThPos]./U[:,:,Param.RhoPos];
  elseif str == "RhoTh"
      @views cOut[:,:,i].=U[:,:,Param.ThPos]
  elseif str == "TPrime"
      p=Pressure(U[:,:,Param.ThPos],U[:,:,Param.ThPos],U[:,:,Param.ThPos],Param);
      @views cOut[:,:,i]=p./(Param.Rd*U[:,:,Param.RhoPos])-Param.TBGrd;
  elseif str == "ThetaPrime"
      @views cOut[:,:,i]=U[:,:,Param.ThPos]./U[:,:,Param.RhoPos]-Param.ThetaBGrd;
  elseif str == "Pres"
      @views cOut[:,:,i]=Pressure(U[:,:,Param.ThPos],U[:,:,Param.ThPos],U[:,:,Param.ThPos],Param);
  elseif str == "Vort"
      @views cOut[:,:,i]=FVort2VecDSS(U[:,:,Param.uPos],U[:,:,Param.vPos],CG,Param);
  end
end
vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
return vtk
end
