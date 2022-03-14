function f = fVelW(x,Param)
switch lower(Param.ProfVelW)
  case 'linear'
    f=x(1);
  otherwise
    f=0;
end
end

