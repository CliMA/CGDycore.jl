function vtkOutput(U,vtkGrid,CG,Param)
nz=Param.Grid.nz;
cOut=zeros(CG.NumG,nz,length(Param.cNames));
for i=1:length(Param.cNames)
  str = Param.cNames[i]
  if str == "Rho"
      cOut[:,:,i]=U[:,:,Param.RhoPos];
  elseif str == "u"
      cOut[:,:,i]=U[:,:,Param.uPos];
  elseif str == "v"
      cOut[:,:,i]=U[:,:,Param.vPos];
  elseif str == "w"
      cOut[:,1,i]=0.5*U[:,1,Param.wPos];
      cOut[:,2:nz-1,i]=0.5*(U[:,1:nz-2,Param.wPos]+U[:,2:nz-1,Param.wPos]);
      cOut[:,nz,i]=0.5*U[:,nz-1,Param.wPos];
  elseif str == "Th"
      cOut[:,:,i]=U[:,:,Param.ThPos]./U[:,:,Param.RhoPos];
  elseif str == "RhoTh"
      cOut[:,:,i]=U[:,:,Param.ThPos]
  elseif str == "TPrime"
      p=Pressure(U[:,:,Param.ThPos],U[:,:,Param.ThPos],U[:,:,Param.ThPos],Param);
      cOut[:,:,i]=p./(Param.Rd*U[:,:,Param.RhoPos])-Param.TBGrd;
  elseif str == "ThetaPrime"
      cOut[:,:,i]=U[:,:,Param.ThPos]./U[:,:,Param.RhoPos]-Param.ThetaBGrd;
  elseif str == "Pres"
      cOut[:,:,i]=Pressure(U[:,:,Param.ThPos],U[:,:,Param.ThPos],U[:,:,Param.ThPos],Param);
  elseif str == "Vort"
      cOut[:,:,i]=FVort2VecDSS(U[:,:,Param.uPos],U[:,:,Param.vPos],CG,Param);
  end
end
vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
return vtk
end
