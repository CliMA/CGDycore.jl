function h=Oro(x,L,H)
p=3;
if x<L/2
  h=H*(x/L)^p;
else
  h=H*((1-x/L))^p;
end
end