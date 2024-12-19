function InitDG2!(Trans,Fcn,Grid,Param)

end

#=
OrdPolyX=Param.OrdPolyX;
nV=Param.nV;
hPos=Param.hPos;
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;


# g=\sum g_i L(\ksi,xQ,i)
# \int L(\ksi,xQ,i) DL(\ksi,xQ,j)abs(DTrans(\ksi)) d\ksi
NX=OrdPolyX+1;
[wQX,ksi]=GaussLobattoQuad(OrdPolyX);
Param.xw=ksi;
Param.DX=DerivativeMatrix(OrdPolyX);
Param.MX=wQX;
Param.DXS=-inv(diag(Param.MX))*Param.DX'*diag(Param.MX);

Param.X=zeros(NX,NX,3,Grid.NumFaces);
Param.dXdx=zeros(NX,NX,3,2,Grid.NumFaces);
Param.dXdxI=zeros(NX,NX,2,3,Grid.NumFaces);
Param.J=zeros(NX,NX,Grid.NumFaces);
Param.rot=zeros(NX,NX,2,3,Grid.NumFaces);
for iF=1:Grid.NumFaces
  for i=1:NX
    for j=1:NX
      [Param.X(i,j,:,iF),Param.J(i,j,iF),Param.dXdx(i,j,:,:,iF)...
        ,Param.dXdxI(i,j,:,:,iF),Param.rot(i,j,:,:,iF)]=...
        Trans(ksi(i),ksi(j),Grid.Faces(iF),Grid);
    end
  end
end

Param.N=zeros(NX,3,Grid.NumEdges);
Param.T1=zeros(NX,3,Grid.NumEdges);
Param.VolSurf=zeros(NX,Grid.NumEdges);
for iE=1:Grid.NumEdges
  iF=Grid.Edges(iE).F(1);
  FE=Grid.Edges(iE).FE(1);
  if FE==1
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(:,1,1,1,iF).^2 ...
                            +Param.dXdx(:,1,2,1,iF).^2 ...
                            +Param.dXdx(:,1,3,1,iF).^2);
    Param.T1(:,:,iE)=Param.dXdx(:,1,:,1,iF)./...
      sqrt(Param.dXdx(:,1,1,1,iF).^2 ...
          +Param.dXdx(:,1,2,1,iF).^2 ....
          +Param.dXdx(:,1,3,1,iF).^2);
    for i=1:NX
      Param.N(i,:,iE)=cross(reshape(Param.T1(i,:,iE),3,1),...
        reshape(Param.X(i,1,:,iF),3,1));
      Param.N(i,:,iE)=Param.N(i,:,iE)/norm(Param.N(i,:,iE),'fro');
    end  
  elseif FE==2
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(NX,:,1,2,iF).^2 ...
                            +Param.dXdx(NX,:,2,2,iF).^2....
                            +Param.dXdx(NX,:,3,2,iF).^2);
    Param.T1(:,:,iE)=Param.dXdx(NX,:,:,2,iF)./...
      sqrt(Param.dXdx(NX,:,1,2,iF).^2 ...
          +Param.dXdx(NX,:,2,2,iF).^2....
          +Param.dXdx(NX,:,3,2,iF).^2);
    for i=1:NX
      Param.N(i,:,iE)=cross(reshape(Param.T1(i,:,iE),3,1),...
        reshape(Param.X(NX,i,:,iF),3,1));
      Param.N(i,:,iE)=Param.N(i,:,iE)/norm(Param.N(i,:,iE),'fro');
    end 
  elseif FE==3
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(:,NX,1,1,iF).^2 ...
                            +Param.dXdx(:,NX,2,1,iF).^2....
                            +Param.dXdx(:,NX,3,1,iF).^2);
    Param.T1(:,:,iE)=Param.dXdx(:,NX,:,1,iF)./...
      sqrt(Param.dXdx(:,NX,1,1,iF).^2 ...
          +Param.dXdx(:,NX,2,1,iF).^2....
          +Param.dXdx(:,NX,3,1,iF).^2);
    for i=1:NX
      Param.N(i,:,iE)=cross(reshape(Param.T1(i,:,iE),3,1),...
        reshape(Param.X(i,NX,:,iF),3,1));
      Param.N(i,:,iE)=Param.N(i,:,iE)/norm(Param.N(i,:,iE),'fro');
    end 
  elseif FE==4
    Param.VolSurf(:,iE)=sqrt(Param.dXdx(1,:,1,2,iF).^2 ...
                            +Param.dXdx(1,:,2,2,iF).^2....
                            +Param.dXdx(1,:,3,2,iF).^2);
    Param.T1(:,:,iE)=Param.dXdx(1,:,:,2,iF)./...
      sqrt(Param.dXdx(1,:,1,2,iF).^2 ...
          +Param.dXdx(1,:,2,2,iF).^2....
          +Param.dXdx(1,:,3,2,iF).^2);
    for i=1:NX
      Param.N(i,:,iE)=cross(reshape(Param.T1(i,:,iE),3,1),...
        reshape(Param.X(1,i,:,iF),3,1));
      Param.N(i,:,iE)=Param.N(i,:,iE)/norm(Param.N(i,:,iE),'fro');
    end
  end
%   Param.N(:,1,iE)=Param.T1(:,2,iE);
%   Param.N(:,2,iE)=-Param.T1(:,1,iE);
end


V=zeros(nV,NX,NX,Grid.NumFaces);
if hPos>0
  V(hPos,:,:,:)=DG2Project(@hF2,Param);
end
if uPos>0
  V(uPos:vPos,:,:,:)=DG2ProjectVec(@VelF2,Param);
  V(uPos,:,:,:)=V(uPos,:,:,:).*V(hPos,:,:,:);
  V(vPos,:,:,:)=V(vPos,:,:,:).*V(hPos,:,:,:);
end
fig=1;
% fig=PlotDG(V(hPos,:,:,:),Grid,Trans,fig,Param);
% fig=PlotDG(V(uPos,:,:,:)./V(hPos,:,:,:),Grid,Trans,fig,Param);
% fig=PlotDG(V(vPos,:,:,:)./V(hPos,:,:,:),Grid,Trans,fig,Param);

time=0;
dtExp=40;
iterEnd=4000;
fig=PlotDG(V(hPos,:,:,:),Grid,Trans,fig,Param);
%fig=PlotDG(V(uPos,:,:,:)./V(hPos,:,:,:),Grid,Trans,fig,Param);
switch Param.IntMethod
  case 'RungeKutta'
    for iter=1:iterEnd
      iter
      V=RungeKuttaExplicit(V,dtExp,Fcn,Grid,Param);
      if mod(iter,1000)==0
        fig=PlotDG(V(hPos,:,:,:),Grid,Trans,fig,Param);
%       fig=PlotDG(V(uPos,:,:,:)./V(hPos,:,:,:),Grid,Trans,fig,Param);
%       fig=PlotDG(V(vPos,:,:,:)./V(hPos,:,:,:),Grid,Trans,fig,Param);
      end
      time=time+dtExp;
    end
end
if Param.OutputEnd
  fig=PlotDG(V(hPos,:,:,:),Grid,Trans,fig,Param);
  fig=PlotDG(V(uPos,:,:,:),Grid,Trans,fig,Param);
end
end
%fig=PlotDG(V(hPos,:,:,:),Grid,Trans,fig,Param);
      
=#

