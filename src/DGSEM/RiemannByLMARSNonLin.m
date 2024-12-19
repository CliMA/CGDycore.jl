function F=RiemannByLMARSNonLin(VLL,VRR,Param)
uPos=Param.uPos;
vPos=Param.vPos;
hPos=Param.hPos;
F=zeros(size(VLL));
pLL=PresShB(VLL,Param);
pRR=PresShB(VRR,Param);
hM=0.5*(VLL(hPos,:)+VRR(hPos,:));
vLL=VLL(uPos,:)./VLL(hPos,:);
vRR=VRR(uPos,:)./VRR(hPos,:);
pM=0.5*(pLL+pRR)...
  -0.5*Param.cS*hM.*(vRR-vLL);
vM = 0.5*(vRR+vLL)...
  -1.0/(2.0*Param.cS).*(pRR-pLL)./hM;
F(hPos,:)= vM.*VLL(hPos,:).*(vM>0)...
            +vM.*VRR(hPos,:).*(vM<0);
          
F(uPos,:)=pM+vM.*VLL(uPos,:).*(vM>0)...
            +vM.*VRR(uPos,:).*(vM<0);
          
F(vPos,:)=   vM.*VLL(vPos,:).*(vM>0)...
            +vM.*VRR(vPos,:).*(vM<0);
end