N = 10;
M = 10;
OrdPoly=3;
H = 400;
hHill=100;
L=1000;
if OrdPoly>0
  %Horizontal Grid
  [wX,xw]=GaussLobattoQuad(OrdPoly);

  Dx=zeros(OrdPoly+1,OrdPoly+1);
  for i=1:OrdPoly+1
    for j=1:OrdPoly+1
      Dx(i,j)=DLagrange(xw(i),xw,j);
    end
  end
  %Vertical Grid
  OrdPolyZ=1;
  [wZ,zw]=GaussLobattoQuad(OrdPolyZ);

  Dz=zeros(OrdPolyZ+1,OrdPolyZ+1);
  for i=1:OrdPolyZ+1
    for j=1:OrdPolyZ+1
      Dz(i,j)=DLagrange(zw(i),zw,j);
    end
  end
end

xP=zeros(M,OrdPoly+1);
xP(1,1)=0;
dx=L/M;
for i=1:M
  for j=1:OrdPoly+1
    xP(i,j)=xP(i,1)+(1+xw(j))/2*dx;
  end
  if i<M
    xP(i+1,1)=xP(i,OrdPoly+1);
  end
end
xP(M,OrdPoly+1)=L;

zP=zeros(M,N+1,OrdPoly+1);
for j=1:M
  for k=1:OrdPoly+1
    zP(j,1,k)=Oro(xP(j,k),L,hHill);
    zP(j,N+1,k)=H;
    dzLoc=(zP(j,N+1,k)-zP(j,1,k))/N;
    for i=2:N
      zP(j,i,k)=zP(j,i-1,k)+dzLoc;
    end
  end
end
zM=zeros(M,N,OrdPoly+1);
for i=1:N
  zM(:,i,:)=0.5*(zP(:,i,:)+zP(:,i+1,:));
end
dz=zeros(M,N,OrdPoly+1);
for i=1:N
  dz(:,i,:)=zP(:,i+1,:)-zP(:,i,:);
end
dzF=zeros(M,N+1,OrdPoly+1);
dzF(:,1,:)=dz(:,1,:);
dzF(:,N+1,:)=dz(:,N,:);
for i=2:N
  dzF(:,i,:)=0.5*(dz(:,i-1,:)+dz(:,i,:));
end
%Metric
ZZ=zeros(OrdPoly+1,2);
X=zeros(M,N,OrdPoly+1,2,2);
J=zeros(M,N,OrdPoly+1,2);
dXdx=zeros(M,N,OrdPoly+1,2,2,2);
dXdxI=zeros(M,N,OrdPoly+1,2,2,2);
dXdxIF=zeros(M,N+1,OrdPoly+1,2,2);
dXdxIC=zeros(M,N,OrdPoly+1,2,2);
JF=zeros(M,N+1,OrdPoly+1);
JC=zeros(M,N,OrdPoly+1);

for i=1:N
  for j=1:M
    ZZ(:,1)=zP(j,i,:);
    ZZ(:,2)=zP(j,i+1,:);
    [X(j,i,:,:,:),J(j,i,:,:),dXdx(j,i,:,:,:,:),dXdxI(j,i,:,:,:,:)]=...
      JacobiDG2(xP(j,:),ZZ,Dx,xw,Dz,zw);
  end
end

for i=1:N
  JC(:,i,:)=0.5*(J(:,i,:,1)+J(:,i,:,2));
  JF(:,i,:)=J(:,i,:,1);
  dXdxIC(:,i,:,:,:)=0.5*(dXdxI(:,i,:,1,:,:)+dXdxI(:,i,:,2,:,:));
  dXdxIF(:,i,:,:,:)=dXdxI(:,i,:,1,:,:);
end
JF(:,N+1,:)=J(:,N,:,2);
dXdxIF(:,N+1,:,:,:)=dXdxI(:,N,:,2,:,:);
u=rand(M,N,OrdPoly+1);
u(:,2:N,:)=0.0;
u=DSS(u,JC);
wF=2.0*rand(M,N+1,OrdPoly+1)-1.0;
wF=DSSF(wF,JC);
wF(:,1,:)=-dXdxIF(:,1,:,2,1).*u(:,1,:)./dXdxIF(:,1,:,2,2);
wF(:,N+1,:)=0.0;
wF=DSSF(wF,JC);

Rho=rand(M,N,OrdPoly+1)+1.0;
RhO(:,:,:)=1;
Rho=DSS(Rho,JC);
RhoF=zeros(M,N+1,OrdPoly+1);
RhoF(:,1,:)=Rho(:,1,:);
RhoF(:,N+1,:)=Rho(:,N,:);
for i=2:N
  RhoF(:,i,:)=(Rho(:,i-1,:).*JC(:,i-1,:)+Rho(:,i,:).*JC(:,i,:))./...
    (JC(:,i-1,:)+JC(:,i,:));
end
K=zeros(M,N,OrdPoly+1);
for i=1:N
  K(:,i,:)=0.5*(u(:,i,:).*u(:,i,:)...
    +0.25*(wF(:,i,:).*wF(:,i,:)+wF(:,i+1,:).*wF(:,i+1,:)));
end
S=zeros(M,N+1,OrdPoly+1);

for i=1:N-1
  uF=(dXdxIC(:,i,:,2,1).*u(:,i,:).*Rho(:,i,:)...
    +dXdxIC(:,i+1,:,2,1).*u(:,i+1,:).*Rho(:,i+1,:))...
    ./(Rho(:,i,:).*dXdxIC(:,i,:,2,2)+Rho(:,i+1,:).*dXdxIC(:,i+1,:,2,2));
  S(:,i+1,:)=RhoF(:,i+1,:).*(wF(:,i+1,:)+uF);
end

%%%%%%%%%%%%%%%%%%
% Part 1
wDotH1 = udwdx(u,wF,Rho,Dx,dXdxIF,JC);
uDotH1 = -wdwdx(wF,Dx,dXdxIF,JC);
uICH1=IntCell(uDotH1.*Rho.*u,JC,wX);
wIFH1=IntFace(wDotH1.*RhoF.*wF,JC,wX);
IH1=uICH1+wIFH1;

%%%%%%%%%%%%%%%%%%
% Part 2
% \rho u \nabla_S K + K \nabla_S \rho u
uDotH2 = dcdx(K,Dx,dXdxIC,JC);
RhoDotH2 = divdx(Rho.*u,Dx,dXdxIC,JC);
uICH2=IntCell(uDotH2.*Rho.*u,JC,wX);
RhoIH2=IntCell(RhoDotH2.*K,JC,wX);
IH2=uICH2+RhoIH2;

%%%%%%%%%%%%%%%%%%%%
% Part 3
wDotV = Sdwdz(wF,RhoF,S,dXdxIF,JC);
uDotV = Sdudz(u,Rho,S,dXdxIF,JC);
RhoDotV = dRhoSdz(S,dXdxIF,JC);
uICV=IntCell(uDotV.*Rho.*u,JC,wX);
wIFV=IntFace(wDotV.*RhoF.*wF,JC,wX);
KICV=IntCell(RhoDotV.*K,JC,wX);
IV=uICV+wIFV+KICV;

%%%%%%%%%%%%%%%%%%%%
% Part 1 + Part 2 + Part 3
uDot=dcdx(K,Dx,dXdxIC,JC)-wdwdx(wF,Dx,dXdxIF,JC)+Sdudz(u,Rho,S,dXdxIF,JC);
wDot=udwdx(u,wF,Rho,Dx,dXdxIF,JC)+Sdwdz(wF,RhoF,S,dXdxIF,JC);
wDot(:,1,:)=0.0;
RhoDot=divdx(Rho.*u,Dx,dXdxIC,JC)+dRhoSdz(S,dXdxIF,JC);
%wDot(:,1,:)=dXdxIF(:,1,:,2,1).*uDot(:,1,:)./dXdxIF(:,1,:,2,2);
uIC=IntCell(uDot.*Rho.*u,JC,wX);
wIF=IntFace(wDot.*RhoF.*wF,JC,wX);
KIC=IntCell(RhoDot.*K,JC,wX);
I=uIC+wIF+KIC;

aa=3;