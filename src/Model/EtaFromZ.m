function eta=EtaFromZ(x,y,z,Param)
etaOld=1.e-7;
F=-Param.Grav*z+PhiBaroWave(x,y,etaOld,Param);
FStart=F;
for i=1:100
  %F=-Param.Grav*z+PhiBaroWave(x,y,etaOld,Param);
  FPrim=-Param.Rd/etaOld*TBaroWave(x,y,etaOld,Param);
  eta=etaOld-F/FPrim;
  if abs(eta-etaOld)<=1.e-13
    break
  end
  etaOld=eta;
  F=-Param.Grav*z+PhiBaroWave(x,y,etaOld,Param);
%   if abs(F)<=abs(FStart)*1.e-12
%     break
%   end
end
end