function StabilityRegion(RK::INT.RungeKuttaMethod;rangeR=-5.0:0.1:2.0,rangeC=-5.0:0.1:5.0)

  Y = zeros(ComplexF64,RK.nStage)
  fY = zeros(ComplexF64,RK.nStage)
  yC = zeros(length(rangeC),length(rangeR))

  for iC in eachindex(rangeC)
    zC = rangeC[iC]  
    for iR in eachindex(rangeR)
      zR = rangeR[iR]  
      z = zR + im * zC
      y0 = complex(1.0)
      @. Y = y0
      for iStage = 1 : RK.nStage
        for jStage = 1 : iStage - 1   
          Y[iStage] += RK.a[iStage,jStage] * fY[jStage]
        end  
        fY[iStage] = z * Y[iStage]  
      end  
      for iStage = 1 : RK.nStage
        y0 += RK.b[iStage] * fY[iStage]
      end  
      yC[iC,iR] = min(abs(y0),1.0)
    end
  end  
  @show size(rangeR),size(rangeC),size(yC)
  p = contour(rangeR,rangeC,yC, levels=10, level_value=0:.1:1 , plot_title=RK.name)
  savefig(p, RK.name*".png")
end

#=
function info=StabilityRegion(Method,NameMethod,mach)
global Param

Param.StenAdvection='Third';
cAMax=3;
cSMax=20;
nA=100;
nS=100;
i=sqrt(-1);
R=zeros(nA+1,nS+1);
nz=40;
switch Method
  case'MIS'
    MIS=MethodMIS(NameMethod);
    for l=1:nz
      z=exp(i*l*pi/(nz+1));
      [RpD,RmD,cAD,cSD]=StabilityMISDiag(nA,nS,cAMax,cSMax,z,MIS);
      contour(cSD,cAD,RpD,[1.00000000001:.05:2.],'linewidth',1.2);
      contour(cSD,cAD,RmD,[1.00000000001:.05:2.],'linewidth',1.2);
%       [Rp,Rm,cA,cS]=StabilityMIS(nA,nS,cAMax,cSMax,z,MIS);
%       contour(cS,cA,Rp,[1.00000000001:.05:2.],'linewidth',1.2);
%       contour(cS,cA,Rm,[1.00000000001:.05:2.],'linewidth',1.2);
      R=max(R,max(RmD,RpD));
    end
  case 'ETD'
    ETD=MethodETD(NameMethod);
    for l=1:20
      z=exp(i*l*pi/21);
      [Rp,Rm,cA,cS]=StabilityETD(nA,nS,cAMax,cSMax,z,ETD);
      contour(cS,cA,Rp,[1.00000001:.05:2.],'linewidth',1.2);
      contour(cS,cA,Rm,[1.00000001:.05:2.],'linewidth',1.2);
      R=max(R,max(Rm,Rp));
    end
end
hold off
contour(cS,cA,R,[1.00000001:.05:2.],'linewidth',1.2);
hold on
plot(cS,mach*cS);
hold off


% RStrang=zeros(nA+1,nS+1);
% Strang=0;
% for l=1:20
%   z=exp(i*l*pi/21);
%   [RStrangp,RStrangm,cA,cS]=StabilityStrang(nA,nS,cAMax,cSMax,z,Strang);
%   contour(cS,cA,RStrangp,[1.00000001:.05:2.],'linewidth',1.2);
%   contour(cS,cA,RStrangm,[1.00000001:.05:2.],'linewidth',1.2);
%   RStrang=max(RStrang,max(RStrangm,RStrangp));
% end
% contour(cS,cA,RStrang,[1.00000001:.05:2.],'linewidth',1.2);
info=0;
end
=#

