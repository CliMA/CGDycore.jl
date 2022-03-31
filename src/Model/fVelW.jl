function fVelW(x,Param)
  str = lower(Param.ProfVelW)
  if str == "linear"
    f=x[1];
  elseif str == "heldsuarezcart"
    f=sin(pi*x[3]/Param.H)
  else
    f=0;
  end
  return f
end

