function f=FluxNonLin(c,iF,Param)
uPos=Param.uPos;
vPos=Param.vPos;
wPos=Param.wPos;
hPos=Param.hPos;
p=PresSh(c(:,:,:),Param);
u=c(uPos,:,:)./c(hPos,:,:);
v=c(vPos,:,:)./c(hPos,:,:);
w=c(wPos,:,:)./c(hPos,:,:);
f=zeros(3,Param.nV,Param.OrdPolyX+1,Param.OrdPolyX+1);
f(1,hPos,:,:)=-c(uPos,:,:);
f(1,uPos,:,:)=-c(uPos,:,:).*u-p;
f(1,vPos,:,:)=-c(uPos,:,:).*v;
f(1,wPos,:,:)=-c(uPos,:,:).*w;

f(2,hPos,:,:)=-c(vPos,:,:);
f(2,uPos,:,:)=-c(vPos,:,:).*u;
f(2,vPos,:,:)=-c(vPos,:,:).*v-p;
f(2,wPos,:,:)=-c(vPos,:,:).*w;


f(3,hPos,:,:)=-c(wPos,:,:);
f(3,uPos,:,:)=-c(wPos,:,:).*u;
f(3,vPos,:,:)=-c(wPos,:,:).*v;
f(3,wPos,:,:)=-c(wPos,:,:).*w-p;
end

