function V=RosenbrockSchur1(V,dt,Fcn,Jac,CG,Param)
Vn=V;
ROS=Param.ROS1;
nV1=size(V,1);
nV2=size(V,2);
nV3=size(V,3);
nJ=nV1*nV2*nV3;
nStage=ROS.nStage;
k=zeros([size(V) nStage]);
[JS]=Jac(V,CG,Param);
if ROS.transformed
  for iStage=1:nStage
    V=Vn;
    for jStage=1:iStage-1
      V=V+ROS.a(iStage,jStage)*k(:,:,:,jStage);
    end
    fV(:,:,:)=Fcn(V,CG,Param);
    for jStage=1:iStage-1
      fV=fV+(ROS.c(iStage,jStage)/dt)*k(:,:,:,jStage);
    end
    k(:,:,:,iStage)=SchurSolve(fV,JS,dt*ROS.d);
  end
  V=Vn;
  for iStage=1:nStage
    V=V+ROS.m(iStage)*k(:,:,:,iStage);
  end
else
  for iStage=1:nStage
    V=Vn;
    for jStage=1:iStage-1
      V=V+ROS.alpha(iStage,jStage)*k(:,:,:,jStage);
    end
    fV(:,:,:)=Fcn(V,CG,Param);
    if iStage>1
      kTemp=zeros(size(Vn));
      for jStage=1:iStage-1
        kTemp=kTemp+ROS.Gamma(iStage,jStage)...
          *k(:,:,:,jStage);
      end
      fV=fV+JacVec(JS,kTemp);
    end
    if ROS.Gamma(iStage,iStage)>0
      fV=fV/ROS.Gamma(iStage,iStage);
      k(:,:,:,iStage)=SchurSolve(fV,JS,dt*ROS.Gamma(iStage,iStage));
    else
      k(:,:,:,iStage)=dt*fV;
    end
  end
  V=Vn;
  for iStage=1:nStage
    V=V+ROS.b(iStage)*k(:,:,:,iStage);
  end
end
end
