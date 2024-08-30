function Interpolation(FaceC,Stencil,Grid,OrdPolynom)
  M = zeros(length(Stencil),size(OrdPolynom,1))
  k = FaceC.Mid
  k /= norm(k)
  @show k
  @show FaceC.Mid
  if abs(k.x) > 0.5 
    t1 = Grids.Point(k.y,-k.x,0.0)
  elseif abs(k.y) > 0.5    
    t1 = Grids.Point(0.0,k.z,-k.y)
  else
    t1.x = -k.z
    t1.y = 0
    t1.z = k.x  
    t1 = Grids.Point(-k.z,0.0,k.x)
  end  
  @show Grids.dot(k,t1)
  t1 /= norm(t1)
  t2 = Grids.cross(t1,k)
  x0 = FaceC.Mid

  for iF in Stencil
    @show iF  
    # Integrate over kites
    for i = 1 :  length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      P2 = Grid.Nodes[iN].P
      if i == 1
        im1 = length(Grid.Faces[iF].N)
      else
        im1 = i - 1
      end  
      iEm1 = Grid.Faces[iF].E[im1]
      P1 = Grid.Edges[iEm1].Mid
      if i == length(Grid.Faces[iF].N)
        ip1 = 1
      else
        ip1 = i + 1
      end
      iEp1 = Grid.Faces[iF].E[ip1]
      P3 = Grid.Edges[iEp1].Mid
      P4 = Grid.Faces[iF].Mid
      @views Integral!(M[iF,:],P1,P2,P3,P4,t1,t2,x0)
    end
  end

#   @. M[i,:] += Integral(P1,P2,P3,P4,
    # Intgrate the monomials
end    

function Integral!(M,P1,P2,P3,P4,t1,t2,x0,OrdPolynom)
  n = 4
  ksi, w = gaussradau(n)
  ksi1 = .5
  ksi2 = .5
  for i = 1 : n
    ksi1 = ksi[i]  
    for j = 1 : n
      ksi2 = ksi[j]  
      P = 0.25 * ((1.0 - ksi1)*(1 - ksi2)*P1 +
                  (1.0 + ksi1)*(1 - ksi2)*P2 +
                  (1.0 + ksi1)*(1 + ksi2)*P3 +
                  (1.0 - ksi1)*(1 + ksi2)*P4)
      A = [P.x -t1.x -t2.x
           P.y -t1.y -t2.y 
           P.z -t1.z -t2.z]
      b = [x0.x     
           x0.y
           x0.z]
      c = A \ b     
      ksi1Hat = c[2]
      ksi2Hat = c[3]
      for k = size(OrdPolynom,1)
        M[k] += w[i] * w[j] * ksi1Hat^OrdPolynom[k,1] * ksi2Hat^OrdPolynom[k,2] 
      end
    end
  end  
end  
