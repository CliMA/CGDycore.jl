function k1F(hL,hC,hR)
  k1F=(hC/(hC+hR))*((hL+hC)/(hL+hC+hR))
end  
function k2F(hL,hC,hR)
  k2F=(hC/(hL+hC))*(hR/(hL+hC+hR))
end

function Rec3(cL,cC,cR,JL,JC,JR)
  kL=(JC/(JC+JL))*((JR+JC)/(JL+JC+JR))
  kR=-(JC/(JR+JC))*(JL/(JL+JC+JR))
  cCL=kL*cL+(1.0-kL-kR)*cC+kR*cR
  kR=(JC/(JC+JR))*((JL+JC)/(JL+JC+JR))
  kL=-(JC/(JL+JC))*(JR/(JL+JC+JR))
  cCR=kL*cL+(1.0-kL-kR)*cC+kR*cR
  return (cCL,cCR)
end  


xLL=0
JL=1 #xLC=.5
xL=1
JC=2 #xC=2 
xR=3
JR=3 #xRC=4.5
xRR=6


#constant function
cL=1.0
cC=1.0
cR=1.0
@show Rec3(cL,cC,cR,JL,JC,JR)

#linear function
cL=.5
cC=2.0
cR=4.5
@show Rec3(cL,cC,cR,JL,JC,JR)

#quadratic function
cL=(xL^3-xLL^3)/JL/3
cC=(xR^3-xL^3)/JC/3
cR=(xRR^3-xR^3)/JR/3
@show cL
@show cC
@show cR
@show Rec3(cL,cC,cR,JL,JC,JR)
@show xL^2
@show xR^2
