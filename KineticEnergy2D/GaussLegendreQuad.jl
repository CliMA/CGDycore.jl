function GaussLegendreQuad(n)
  w=zeros(n)
  x=zeros(n)
  if n == 1
    w[1]=2
    x[1]=0
  elseif n == 2
    w[1]=1
    w[2]=1
    x[1]=-1/sqrt(3)
    x[2]=1/sqrt(3)
  elseif n == 3
    w[2]=8/9
    w[1]=5/9
    w[3]=5/9
    
    x[2]=0
    x[1]=-sqrt(3/5)
    x[3]=sqrt(3/5)
    
  elseif n == 4
    w[2]=(18+sqrt(30))/36
    w[3]=(18+sqrt(30))/36
    w[1]=(18-sqrt(30))/36
    w[4]=(18-sqrt(30))/36
    
    x[2]=-sqrt(3/7-2/7*sqrt(6/5))
    x[3]=sqrt(3/7-2/7*sqrt(6/5))
    x[1]=-sqrt(3/7+2/7*sqrt(6/5))
    x[4]=sqrt(3/7+2/7*sqrt(6/5))
    
  elseif n == 5
    w[3]=128/225
    w[2]=(322+13*sqrt(70))/900
    w[4]=(322+13*sqrt(70))/900
    w[1]=(322-13*sqrt(70))/900
    w[5]=(322-13*sqrt(70))/900
    
    x[3]=0
    x[2]=-1/3*sqrt(5-2*sqrt(10/7))
    x[4]=1/3*sqrt(5-2*sqrt(10/7))
    x[1]=-1/3*sqrt(5+2*sqrt(10/7))
    x[5]=1/3*sqrt(5+2*sqrt(10/7))
  else
    println("ord1+ord2 too large or small")
  end
  return (w,x)
end

