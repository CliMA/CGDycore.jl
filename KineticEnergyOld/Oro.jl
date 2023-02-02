function Oro(x,L,H)
  p=3
  if x<L/2
    h=H*(2*x/L)^p
  else
    h=H*(2*(1-x/L))^p
  end
  return h
end
