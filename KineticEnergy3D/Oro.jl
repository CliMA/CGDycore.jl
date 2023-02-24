function Oro(x,y,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    h = 0.0  
  elseif Example == "HillAgnesiCart"
    h = Param.h / (1.0 + ((x - Param.xc)/Param.a)^2)  
  else  
    p=3
    if x<L/2
      h=H*(2*x/L)^p
    else
      h=H*(2*(1-x/L))^p
    end  
  end
  return h
end
