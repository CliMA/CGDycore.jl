function Oro(x,y,Param)
  Example = Param.Example
  if Example == "WarmBubble2DXCart"
    h = 0.0  
    h = 4000 / (1.0 + ((x - 10000.0)/2000.0)^2)  
  elseif Example == "HillAgnesiXCart" 
    h = Param.h / (1.0 + ((x - Param.xc)/Param.a)^2)  
  elseif Example == "HillAgnesiYCart"
    h = Param.h / (1.0 + ((y - Param.yc)/Param.a)^2)  
  elseif Example == "AdvectionCart"
    h = 0.0
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
