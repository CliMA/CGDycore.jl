function ProjectW(Fun,time,CG,Metric,Global,Param)
OrdPoly=CG.OrdPoly;
nz=Global.Grid.nz;
w=zeros(nz,CG.NumG);
fLoc=zeros(OrdPoly+1,OrdPoly+1);
X = Metric.X
JF = Metric.JF
for iz=1:nz-1
  for iF=1:Global.Grid.NumFaces
    for j=1:OrdPoly+1
      for i=1:OrdPoly+1
        x=0.5*(X[i,j,2,:,iz,iF]+X[i,j,1,:,iz+1,iF]);
        w[iz,CG.Glob[i,j,iF]]=Fun(x,time,Global,Param)
      end
    end
  end
end
return w
end

function ProjectW!(w,Fun,time,CG,Metric,Global,Param)
  OrdPoly=CG.OrdPoly;
  nz=Global.Grid.nz;
  X = Metric.X
  JF = Metric.JF
  @. w = 0.0
  @inbounds for iz=1:nz-1
    @inbounds for iF=1:Global.Grid.NumFaces
      @inbounds for j=1:OrdPoly+1
        @inbounds for i=1:OrdPoly+1
          x=0.5*(X[i,j,2,:,iz,iF]+X[i,j,1,:,iz+1,iF])
          ind = CG.Glob[i,j,iF]
          @show x
          w[iz,ind]=Fun(x,time,Global,Param)
        end
      end
    end
  end
end
