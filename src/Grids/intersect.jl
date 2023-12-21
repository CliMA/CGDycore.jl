function intersect(Face1, Grid1, Face2, Grid2, pointsSameEPS)

  intersector = Face1
  intersection = Face2
  intor_id = "P1"
  inted_id = "P2"
  intersectorFace1 = Face1.Area > Face2.Area
  if isapprox(Face1.Area, Face2.Area)
    dc = Face1.Mid - Face2.Mid
    if dc.x > 0. || (dc.x == 0. && (dc.y > 0. || (dc.y == 0. && dc.z > 0.)))
      intersectorFace1 = true
    end
  end
  if intersectorFace1
    intersector = Polygon(Face1,Grid1.Nodes)
    intersection = Polygon(Face2,Grid2.Nodes)
  else 
    intor_id = "P2"
    inted_id = "P1"
    intersector = Polygon(Face2,Grid2.Nodes)
    intersection = Polygon(Face1,Grid1.Nodes)
  end  
  if intersection.valid
    for i = 1 : intersector.size
      s1 = intersector.P[i]
      if i < intersector.size
        s2 = intersector.P[i+1] 
      else  
        s2 = intersector.P[1]
      end  
      clip!(intersection,GreatCircle(s1,s2),pointsSameEPS)
    end
  end
  if intersection.size > 2
    triangles = triangulate(intersection)
    A = area(triangles)
  else
    A = 0.0
  end
  return A
end  

function clip!(Pg,C,pointsSameEPS)
  ClippedPg = Polygon()
  first_in = inLeftHemisphere(Pg.P[1],C,-1.5*EPS)

  for i = 1 : Pg.size 
    in = i % Pg.size + 1
    second_in = inLeftHemisphere(Pg.P[in],C,-1.5*EPS)
    if first_in && second_in
      ClippedPg.size += 1
      ClippedPg.P[ClippedPg.size] = Pg.P[in]  
    elseif ~first_in && ~second_in
      continue
    else
      segment = GreatCircle(Pg.P[i],Pg.P[in])  
      ip = intersect(C,segment,pointsSameEPS)
      if ip.x == 1 && ip.y == 1 && ip.z == 1
        # consider the segments parallel
        ClippedPg.size += 1
        ClippedPg.P[ClippedPg.size] = Pg.P[in]  
        first_in = second_in
        continue
      end
      if second_in
        inn = (in+1) % Pg.size + 1
        segment_n = GreatCircle(Pg.P[in], Pg.P[inn])
        if inLeftHemisphere(ip,segment,-1.5*EPS) &&
          inLeftHemisphere(ip,segment_n,-1.5*EPS) &&
          distance(ip,Pg.P[in]) > pointsSameEPS
          ClippedPg.size += 1
          ClippedPg.P[ClippedPg.size] = ip
        end  
        ClippedPg.size += 1
        ClippedPg.P[ClippedPg.size] =Pg.P[in]  
      else  
        if distance(ip, Pg.P[i]) > pointsSameEPS
          ClippedPg.size += 1
          ClippedPg.P[ClippedPg.size] = ip
        end    
      end    
    end    
    first_in = second_in
  end
  # Update polygon
  Pg.size = ClippedPg.size
  if Pg.size < 3
    Pg.valid = false  
  else  
    for i = 1 : Pg.size
      Pg.P[i] = ClippedPg.P[i]
    end
  end  
end

function intersect(C1::GreatCircle,C2::GreatCircle,pointsSameEPS)
  sp = cross(cross(C1.P1,C1.P2),cross(C2.P1,C2.P2))
  sp_norm = norm(sp)
  gcircles_distinct = sp_norm > EPS
  if gcircles_distinct
    sp /= sp_norm  
    sp_2 = -sp
    d1 = max(distance2(sp, C2.P1), distance2(sp, C2.P2));
    d2 = max(distance2(sp_2, C2.P1), distance2(sp_2, C2.P2));
    if d1 < d2
      return sp  
    else
      return sp_2  
    end
  else
    return Point(1,1,1)
  end  
end




