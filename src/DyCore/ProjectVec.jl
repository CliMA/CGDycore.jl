function ProjectVec(Fun,time,CG,Metric,Phys,Global,Param)
  OrdPoly=CG.OrdPoly;
  OrdPolyZ=CG.OrdPolyZ;
  nz=Global.Grid.nz;
  uS=zeros(nz,CG.NumG);
  vS=zeros(nz,CG.NumG);
  uLoc=zeros(OrdPoly+1,OrdPoly+1);
  vLoc=zeros(OrdPoly+1,OrdPoly+1);
  X = Metric.X
  for iF=1:Global.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        ind = CG.Glob[i,j,iF]
        for iz=1:nz
          x=0.5*(X[i,j,1,:,iz,iF]+X[i,j,2,:,iz,iF]);
          (uLoc,vLoc)=Fun(x,time,Phys,Global,Param);
          uS[iz,ind]=uLoc
          vS[iz,ind]=vLoc
        end
      end
    end
  end
  return (uS,vS)
end

function ProjectVec!(uS,vS,Fun,time,CG,Metric,Phys,Global,Param)
  OrdPoly=CG.OrdPoly;
  OrdPolyZ=CG.OrdPolyZ;
  nz=Global.Grid.nz;
  X = Metric.X
  @. uS = 0.0
  @. vS = 0.0
  x=zeros(3)
  uLoc = 0.0
  vLoc = 0.0
  @inbounds for iF=1:Global.Grid.NumFaces
    @inbounds for j=1:OrdPoly+1
      @inbounds for i=1:OrdPoly+1
        ind = CG.Glob[i,j,iF]
        @inbounds for iz=1:nz
          @views @. x = 0.5*(X[i,j,1,:,iz,iF]+X[i,j,2,:,iz,iF]);
          (uLoc,vLoc)=Fun(x,time,Phys,Global,Param);
          uS[iz,ind] = uLoc
          vS[iz,ind] = vLoc
        end
      end
    end
  end
end

