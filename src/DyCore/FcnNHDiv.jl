function [divc,divuF,divwF]=FDivCompEuler(c,cCG,uF,uFCG,wF,wFCG,CG,Param)
OrdPolyX=CG.OrdPolyX;
OrdPolyY=CG.OrdPolyY;
divc=zeros(size(c));
divuF=zeros(size(uF));
divwF=zeros(size(wF));
nV=4;
uLoc=zeros(OrdPolyX+1,OrdPolyY+1,Param.Grid.NumFaces,4);
divLoc=zeros(OrdPolyX+1,OrdPolyY+1,Param.Grid.NumFaces,4);
RhoPos=1;
uPos=2;
wPos=3;
ThPos=4;

for iF=1:Param.Grid.NumFaces
  uLoc(:,:,iF,uPos)=uFCG.IntXLGL*uF(:,:,iF)*uFCG.IntYLGL';
  uLoc(:,:,iF,wPos)=wFCG.IntXLGL*wF(:,:,iF)*wFCG.IntYLGL';
  uLoc(:,:,iF,ThPos)=cCG.IntXLGL*c(:,:,iF,Param.ThPos)*cCG.IntYLGL';
  uLoc(:,:,iF,RhoPos)=cCG.IntXLGL*c(:,:,iF,Param.RhoPos)*cCG.IntYLGL';
  pLoc=Pressure(uLoc(:,:,iF,ThPos),Param);
 
  for iV=1:nV
  divLoc(:,:,iF,iV)=-CG.DWX*((-reshape(CG.dXdx(:,:,1,2,iF),OrdPolyX+1,OrdPolyY+1).*uLoc(:,:,iF,wPos)+...
    reshape(CG.dXdx(:,:,2,2,iF),OrdPolyX+1,OrdPolyY+1).*uLoc(:,:,iF,uPos))...
    .*(uLoc(:,:,iF,iV)./uLoc(:,:,iF,RhoPos)))...
    -((-reshape(CG.dXdx(:,:,2,1,iF),OrdPolyX+1,OrdPolyY+1).*uLoc(:,:,iF,uPos)+...
    reshape(CG.dXdx(:,:,1,1,iF),OrdPolyX+1,OrdPolyY+1).*uLoc(:,:,iF,wPos))...
    .*(uLoc(:,:,iF,iV)./uLoc(:,:,iF,RhoPos)))*CG.DWY';
  end
  divLoc(:,:,iF,uPos)=divLoc(:,:,iF,uPos)-...
    reshape(CG.dXdx(:,:,2,2,iF),OrdPolyX+1,OrdPolyY+1).*(CG.DSX*pLoc)+...
    reshape(CG.dXdx(:,:,2,1,iF),OrdPolyX+1,OrdPolyY+1).*(pLoc*CG.DSY');
  divLoc(:,:,iF,wPos)=divLoc(:,:,iF,wPos)+...
    reshape(CG.dXdx(:,:,1,2,iF),OrdPolyX+1,OrdPolyY+1).*(CG.DSX*pLoc)-...
    reshape(CG.dXdx(:,:,1,1,iF),OrdPolyX+1,OrdPolyY+1).*(pLoc*CG.DSY');
  
  
  if Param.Coriolis
    divLoc(:,:,iF,uPos)=divLoc(:,:,iF,uPos)+2.0*Param.Omega*sin(CG.lat(:,:,iF)).*CG.J(:,:,iF).*uLoc(:,:,iF,wPos);
    divLoc(:,:,iF,wPos)=divLoc(:,:,iF,wPos)-2.0*Param.Omega*sin(CG.lat(:,:,iF)).*CG.J(:,:,iF).*uLoc(:,:,iF,uPos);
  end
  if Param.Gravitation
    divLoc(:,:,iF,wPos)=divLoc(:,:,iF,wPos)-Param.Grav*CG.J(:,:,iF).*uLoc(:,:,iF,RhoPos);
  end
end

for iE=1:Param.Grid.NumEdges
  if strcmp(Param.Grid.Edges(iE).Type,'X')
    if size(Param.Grid.Edges(iE).F,2)==2
      iF1=Param.Grid.Edges(iE).F(1);
      SE1=Param.Grid.Edges(iE).FE(1);
      iF2=Param.Grid.Edges(iE).F(2);
      SE2=Param.Grid.Edges(iE).FE(2);
      LenEdge=CG.OrdPolyX+1;
      VLL=reshape(uLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,:),LenEdge,nV);
      vTemp=CG.N(iE).N(:,1).*VLL(:,uPos)+CG.N(iE).N(:,2).*VLL(:,wPos);
      VLL(:,wPos)=CG.T1(iE).T(:,1).*VLL(:,uPos)+CG.T1(iE).T(:,2).*VLL(:,wPos);
      VLL(:,uPos)=vTemp;
      VRR=reshape(uLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,:),LenEdge,nV);
      vTemp=CG.N(iE).N(:,1).*VRR(:,uPos)+CG.N(iE).N(:,2).*VRR(:,wPos);
      VRR(:,wPos)=CG.T1(iE).T(:,1).*VRR(:,uPos)+CG.T1(iE).T(:,2).*VRR(:,wPos);
      VRR(:,uPos)=vTemp;
      FLoc=RiemannL(VLL,VRR,Param);
      Temp=FLoc(:,uPos);
      FLoc(:,uPos)=CG.N(iE).N(:,1).*Temp+CG.T1(iE).T(:,1).*FLoc(:,wPos);
      FLoc(:,wPos)=CG.N(iE).N(:,2).*Temp+CG.T1(iE).T(:,2).*FLoc(:,wPos);
      for iV=1:nV
        temp=FLoc(:,iV).*CG.VolSurf(iE).VolSurf;
        divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
          =divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
          -reshape(temp,length(CG.i1(SE1).Ind),length(CG.i2(SE1).Ind))/CG.wY(1);
      end
      
      FLoc=RiemannR(VLL,VRR,Param);
      Temp=FLoc(:,uPos);
      FLoc(:,uPos)=CG.N(iE).N(:,1).*Temp+CG.T1(iE).T(:,1).*FLoc(:,wPos);
      FLoc(:,wPos)=CG.N(iE).N(:,2).*Temp+CG.T1(iE).T(:,2).*FLoc(:,wPos);
      for iV=1:nV
        temp=FLoc(:,iV).*CG.VolSurf(iE).VolSurf;
        divLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,iV)...
          =divLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,iV)...
          +reshape(temp,length(CG.i1(SE2).Ind),length(CG.i2(SE2).Ind))/CG.wY(end);
      end
    else
      %Boundary
      iF1=Param.Grid.Edges(iE).F(1);
      SE1=Param.Grid.Edges(iE).FE(1);
      VLL=reshape(uLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,:),OrdPolyX+1,nV);
      vTemp=CG.N(iE).N(:,1).*VLL(:,uPos)+CG.N(iE).N(:,2).*VLL(:,wPos);
      VLL(:,wPos)=CG.T1(iE).T(:,1).*VLL(:,uPos)+CG.T1(iE).T(:,2).*VLL(:,wPos);
      VLL(:,uPos)=vTemp;
      VRR=VLL;
      if SE1==1 || SE1==2
        VRR(:,uPos)=-VLL(:,uPos);
      else
        VRR(:,uPos)=VLL(:,uPos);
        VLL(:,uPos)=-VRR(:,uPos);
      end
      VRR(:,wPos)=VLL(:,wPos);
      FLoc=RiemannL(VLL,VRR,Param);
      Temp=FLoc(:,uPos);
      FLoc(:,uPos)=CG.N(iE).N(:,1).*Temp+CG.T1(iE).T(:,1).*FLoc(:,wPos);
      FLoc(:,wPos)=CG.N(iE).N(:,2).*Temp+CG.T1(iE).T(:,2).*FLoc(:,wPos);
      if SE1==1
        for iV=1:nV
          FLoc(:,iV)=FLoc(:,iV).*CG.VolSurf(iE).VolSurf;
          divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
            =divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
            -FLoc(:,iV)/CG.wY(1);
        end
      else
        for iV=1:nV
          FLoc(:,iV)=FLoc(:,iV).*CG.VolSurf(iE).VolSurf;
          divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
            =divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
            +FLoc(:,iV)/CG.wY(end);
        end
      end
    end
  elseif strcmp(Param.Grid.Edges(iE).Type,'Y')   
    if size(Param.Grid.Edges(iE).F,2)==2
      iF1=Param.Grid.Edges(iE).F(1);
      SE1=Param.Grid.Edges(iE).FE(1);
      iF2=Param.Grid.Edges(iE).F(2);
      SE2=Param.Grid.Edges(iE).FE(2);
      VLL=reshape(uLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,:),OrdPolyY+1,nV);
      vTemp=CG.N(iE).N(:,1).*VLL(:,uPos)+CG.N(iE).N(:,2).*VLL(:,wPos);
      VLL(:,wPos)=CG.T1(iE).T(:,1).*VLL(:,uPos)+CG.T1(iE).T(:,2).*VLL(:,wPos);
      VLL(:,uPos)=vTemp;
      VRR=reshape(uLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,:),OrdPolyY+1,nV);
      vTemp=CG.N(iE).N(:,1).*VRR(:,uPos)+CG.N(iE).N(:,2).*VRR(:,wPos);
      VRR(:,wPos)=CG.T1(iE).T(:,1).*VRR(:,uPos)+CG.T1(iE).T(:,2).*VRR(:,wPos);
      VRR(:,uPos)=vTemp;
      FLoc=RiemannL(VLL,VRR,Param);
      Temp=FLoc(:,uPos);
      FLoc(:,uPos)=CG.N(iE).N(:,1).*Temp+CG.T1(iE).T(:,1).*FLoc(:,wPos);
      FLoc(:,wPos)=CG.N(iE).N(:,2).*Temp+CG.T1(iE).T(:,2).*FLoc(:,wPos);
      for iV=1:nV
        FLoc(:,iV)=FLoc(:,iV).*CG.VolSurf(iE).VolSurf;
        divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
          =divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
          -FLoc(:,iV)'/CG.wX(1);
%         divLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,iV)...
%           =divLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,iV)...
%           +FLoc(:,iV)'/CG.wX(end);
      end
      FLoc=RiemannR(VLL,VRR,Param);
      Temp=FLoc(:,uPos);
      FLoc(:,uPos)=CG.N(iE).N(:,1).*Temp+CG.T1(iE).T(:,1).*FLoc(:,wPos);
      FLoc(:,wPos)=CG.N(iE).N(:,2).*Temp+CG.T1(iE).T(:,2).*FLoc(:,wPos);
      for iV=1:nV
        FLoc(:,iV)=FLoc(:,iV).*CG.VolSurf(iE).VolSurf;
%         divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
%           =divLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,iV)...
%           -FLoc(:,iV)'/CG.wX(1);
        divLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,iV)...
          =divLoc(CG.i1(SE2).Ind,CG.i2(SE2).Ind,iF2,iV)...
          +FLoc(:,iV)'/CG.wX(end);
      end
    else
      %Boundary
      iF1=Param.Grid.Edges(iE).F(1);
      SE1=Param.Grid.Edges(iE).FE(1);
      VLL=reshape(uLoc(CG.i1(SE1).Ind,CG.i2(SE1).Ind,iF1,:),OrdPolyY+1,nV);
      vTemp=CG.N(iE).N(:,1).*VLL(:,uPos)+CG.N(iE).N(:,2).*VLL(:,wPos);
      VLL(:,wPos)=CG.T1(iE).T(:,1).*VLL(:,uPos)+CG.T1(iE).T(:,2).*VLL(:,wPos);
      VLL(:,uPos)=vTemp;
      VRR=VLL;
      if SE1==1 || SE1==2
        VRR(:,uPos)=-VLL(:,uPos);
      else
        VRR(:,uPos)=VLL(:,uPos);
        VLL(:,uPos)=-VRR(:,uPos);
      end
      VRR(:,wPos)=VLL(:,wPos);
    end
    
  end
end
for iF=1:Param.Grid.NumFaces
  divc(:,:,iF,Param.ThPos)=cCG.IntXLG*divLoc(:,:,iF,ThPos)*cCG.IntYLG';
  divc(:,:,iF,Param.RhoPos)=cCG.IntXLG*divLoc(:,:,iF,RhoPos)*cCG.IntYLG';
  divuF(:,:,iF,Param.uPos)=uFCG.IntXLG*divLoc(:,:,iF,uPos)*uFCG.IntYLG';
  divwF(:,:,iF,Param.wPos)=wFCG.IntXLG*divLoc(:,:,iF,wPos)*wFCG.IntYLG';
end
end


function F=Riemann(VL,VR,Param)
pL=Pressure(VL(:,4),Param);
pR=Pressure(VR(:,4),Param);
vL=VL(:,2)./VL(:,1);
vR=VR(:,2)./VR(:,1);
wL=VL(:,3)./VL(:,1);
wR=VR(:,3)./VR(:,1);
thL=VL(:,4)./VL(:,1);
thR=VR(:,4)./VR(:,1);

F=zeros(size(VL));
VM=0.5*(VL(:,2)+VR(:,2));
pM=0.5*(pL+pR);
F(:,1)=VM;
for i=1:size(VL,1)
  if VM(i)>0
    F(i,2)=pM(i)+VM(i).*vL(i);
    F(i,3)=VM(i).*wL(i);
    F(i,4)=VM(i).*thL(i);
  else
    F(i,2)=pM(i)+VM(i)*vR(i);
    F(i,3)=VM(i).*wR(i);
    F(i,4)=VM(i).*thR(i);
  end
end
end

function F=RiemannL(VL,VR,Param)
pL=Pressure(VL(:,4),Param);
pR=Pressure(VR(:,4),Param);
vL=VL(:,2)./VL(:,1);
vR=VR(:,2)./VR(:,1);
wL=VL(:,3)./VL(:,1);
wR=VR(:,3)./VR(:,1);
thL=VL(:,4)./VL(:,1);
thR=VR(:,4)./VR(:,1);

F=zeros(size(VL));
VM=0.5*(VL(:,2)+VR(:,2));
pM=0.5*(pL+pR);
F(:,1)=VM;
for i=1:size(VL,1)
  if VM(i)>0
    F(i,2)=(pM(i)-pL(i))+VM(i).*vL(i);
    F(i,3)=VM(i).*wL(i);
    F(i,4)=VM(i).*thL(i);
  else
    F(i,2)=(pM(i)-pL(i))+VM(i)*vR(i);
    F(i,3)=VM(i).*wR(i);
    F(i,4)=VM(i).*thR(i);
  end
end
end

function F=RiemannR(VL,VR,Param)
pL=Pressure(VL(:,4),Param);
pR=Pressure(VR(:,4),Param);
vL=VL(:,2)./VL(:,1);
vR=VR(:,2)./VR(:,1);
wL=VL(:,3)./VL(:,1);
wR=VR(:,3)./VR(:,1);
thL=VL(:,4)./VL(:,1);
thR=VR(:,4)./VR(:,1);

F=zeros(size(VL));
VM=0.5*(VL(:,2)+VR(:,2));
pM=0.5*(pL+pR);
F(:,1)=VM;
for i=1:size(VL,1)
  if VM(i)>0
    F(i,2)=(pM(i)-pR(i))+VM(i).*vL(i);
    F(i,3)=VM(i).*wL(i);
    F(i,4)=VM(i).*thL(i);
  else
    F(i,2)=(pM(i)-pR(i))+VM(i)*vR(i);
    F(i,3)=VM(i).*wR(i);
    F(i,4)=VM(i).*thR(i);
  end
end
end


