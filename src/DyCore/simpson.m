function res=simpson(x0,xN,r,dx,f)
  n=(xN-x0)/dx+1;
  h=(xN-x0)/n;
  res=0.5*(f(x0,r)+f(xN,r));
  xi=x0;
  for i=1:(n-1)
    xi=xi+h;
    res=res+f(xi,r);
  end
  xi=x0-0.5*h;
  for i=1:n
    xi=xi+h;
    res=res+2.0*f(xi,r);
  end
  res=res*h/3.0;
end