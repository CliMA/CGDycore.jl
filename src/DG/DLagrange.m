function Df=DLagrange(x,xP,iP)
Df=zeros(size(x));
f=ones(size(x));
for i=1:size(xP,1)
   if i~=iP
      Df=Df.*(x-xP(i))./(xP(iP)-xP(i))+f./(xP(iP)-xP(i));
      f=f.*(x-xP(i))./(xP(iP)-xP(i));
   end
end
end