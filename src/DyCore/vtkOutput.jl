function vtk=vtkOutput(U,vtkGrid,CG,Param)
nz=Param.Grid.nz;
cOut=zeros(CG.NumG,nz,size(Param.cNames,2));
for i=1:size(Param.cNames,2)
  switch Param.cNames(i).s
    case 'Rho'
      cOut(:,:,i)=U(:,:,Param.RhoPos);
    case 'u'
      cOut(:,:,i)=U(:,:,Param.uPos);
    case 'w'
      cOut(:,1,i)=0.5*U(:,1,Param.wPos);
      cOut(:,2:nz-1,i)=0.5*(U(:,1:nz-2,Param.wPos)+U(:,2:nz-1,Param.wPos));
      cOut(:,nz,i)=0.5*U(:,nz-1,Param.wPos);
    case 'Th'
      cOut(:,:,i)=U(:,:,Param.ThPos)./U(:,:,Param.RhoPos);
    case 'TPrime'
      p=Pressure(U(:,:,Param.ThPos),U(:,:,Param.ThPos),U(:,:,Param.ThPos),Param);
      cOut(:,:,i)=p./(Param.Rd*U(:,:,Param.RhoPos))-Param.TBGrd;
    case 'ThetaPrime'
      cOut(:,:,i)=U(:,:,Param.ThPos)./U(:,:,Param.RhoPos)-Param.ThetaBGrd;
    case 'Pres'
      cOut(:,:,i)=Pressure(U(:,:,Param.ThPos),U(:,:,Param.ThPos),U(:,:,Param.ThPos),Param);
  end
end
vtk=vtkCG(cOut,CG,Param,vtkGrid,Param.vtk);
end
