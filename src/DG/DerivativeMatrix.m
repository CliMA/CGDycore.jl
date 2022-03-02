function [CG]=DerivativeMatrix(CG)
OrdPolyX=CG.OrdPolyX;
CG.DSX=zeros(OrdPolyX+1,OrdPolyX+1);
for i=1:OrdPolyX+1
  for j=1:OrdPolyX+1
    CG.DSX(i,j)=DLagrange(CG.xwX(i),CG.xwX,j);
  end
end
CG.DWX=-inv(diag(CG.wX))*CG.DSX'*diag(CG.wX);

OrdPolyY=CG.OrdPolyY;
CG.DSY=zeros(OrdPolyY+1,OrdPolyY+1);
for i=1:OrdPolyY+1
  for j=1:OrdPolyY+1
    CG.DSY(i,j)=DLagrange(CG.xwY(i),CG.xwY,j);
  end
end
CG.DWY=-inv(diag(CG.wY))*CG.DSY'*diag(CG.wY);

if strcmp(CG.FeX,'LG')
  OrdPolyX=CG.OrdPolyX;
  CG.DSXLGL=zeros(OrdPolyX+2,OrdPolyX+1);
  CG.IntXLGL=zeros(OrdPolyX+2,OrdPolyX+1);
  CG.IntXLG=zeros(OrdPolyX+1,OrdPolyX+2);
  [CG.wXLGL,CG.xwXLGL]=GaussLobattoQuad(CG.OrdPolyX+1);
  for i=1:OrdPolyX+2
    for j=1:OrdPolyX+1
      CG.DSXLGL(i,j)=DLagrange(CG.xwXLGL(i),CG.xwX,j);
      CG.IntXLGL(i,j)=Lagrange(CG.xwXLGL(i),CG.xwX,j);
      CG.IntXLG(j,i)=Lagrange(CG.xwX(j),CG.xwXLGL,i);
    end
  end
  %CG.IntXLG=inv(diag(CG.wX))*CG.IntXLG*diag(CG.wXLGL);
  CG.DWXLGL=-inv(diag(CG.wX))*CG.DSXLGL'*diag(CG.wXLGL);
end
if strcmp(CG.FeX,'LG') && strcmp(CG.FeY,'LGL')
  DSXLGL=zeros(OrdPolyX+2,OrdPolyX+1);
  for i=1:OrdPolyX+2
    for j=1:OrdPolyX+1
      DSXLGL(i,j)=DLagrange(CG.xwXLGL(i),CG.xwX,j);
    end
  end
  CG.DWXLG=-inv(diag(CG.wX))*DSXLGL'*diag(CG.wXLGL);
end
if strcmp(CG.FeX,'LGL')
  [CG.wXLG,CG.xwXLG]=GaussLegendreQuad(CG.OrdPolyX-1);
  CG.IntXLG=zeros(CG.OrdPolyX,CG.OrdPolyX+1);
  for i=1:CG.OrdPolyX
    for j=1:CG.OrdPolyX+1
      CG.IntXLG(i,j)=Lagrange(CG.xwXLG(i),CG.xwX,j);
      CG.IntXLGL(j,i)=Lagrange(CG.xwX(j),CG.xwXLG,i);
    end
  end
  CG.IntXLGL=eye(CG.OrdPolyX+1);
  CG.IntXLG=eye(CG.OrdPolyX+1);
end

if strcmp(CG.FeY,'LGL')
  [CG.wYLG,CG.xwYLG]=GaussLegendreQuad(CG.OrdPolyY-1);
  CG.IntYLG=zeros(CG.OrdPolyY,CG.OrdPolyY+1);
  for i=1:CG.OrdPolyY
    for j=1:CG.OrdPolyY+1
      CG.IntYLG(i,j)=Lagrange(CG.xwYLG(i),CG.xwY,j);
      CG.IntYLGL(j,i)=Lagrange(CG.xwY(j),CG.xwYLG,i);
    end
  end
  CG.IntYLGL=eye(CG.OrdPolyY+1);
  CG.IntYLG=eye(CG.OrdPolyY+1);
end

if strcmp(CG.FeY,'LG')
  OrdPolyY=CG.OrdPolyY;
  CG.DSYLGL=zeros(OrdPolyY+2,OrdPolyY+1);
  CG.IntYLGL=zeros(OrdPolyY+2,OrdPolyY+1);
  CG.IntYLG=zeros(OrdPolyY+1,OrdPolyY+2);
  [CG.wYLGL,CG.xwYLGL]=GaussLobattoQuad(CG.OrdPolyY+1);
  for i=1:OrdPolyY+2
    for j=1:OrdPolyY+1
      CG.DSYLGL(i,j)=DLagrange(CG.xwYLGL(i),CG.xwY,j);
      CG.IntYLGL(i,j)=Lagrange(CG.xwYLGL(i),CG.xwY,j);
      CG.IntYLG(j,i)=Lagrange(CG.xwY(j),CG.xwYLGL,i);
    end
  end
  %CG.IntYLG=inv(diag(CG.wY))*CG.IntYLG*diag(CG.wYLGL);
  CG.DWYLGL=-inv(diag(CG.wY))*CG.DSYLGL'*diag(CG.wYLGL);
end
end
