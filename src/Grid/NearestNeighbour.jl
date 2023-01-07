function InsideFace(P,Face,Grid)

# Points in polygon oriented, such that the normal is directed outside

# (P1-M)x(P2-M)*(P-M)
  P1 = Grid.Nodes[Face.N[1]].P
  P2 = Grid.Nodes[Face.N[2]].P
  det1 = dot(cross(P1,P2),P)
  P1 = Grid.Nodes[Face.N[2]].P
  P2 = Grid.Nodes[Face.N[3]].P
  det2 = dot(cross(P1,P2),P)
  P1 = Grid.Nodes[Face.N[3]].P
  P2 = Grid.Nodes[Face.N[4]].P
  det3 = dot(cross(P1,P2),P)
  P1 = Grid.Nodes[Face.N[4]].P
  P2 = Grid.Nodes[Face.N[1]].P
  det4 = dot(cross(P1,P2),P)
  if det1 >= 0.0 && det2 >= 0.0 && det3 >= 0.0 && det4 >= 0.0
    return true  
  else
    return false
  end  

end

function walk_to_nc(target_cc,start_Face,xw,Trans,Rad,Grid)

  Face_cc =  Grid.Faces[start_Face].Mid
  searching = true
  # calculate a measure for the distance to target point
  sp = dot(target_cc, Face_cc)
  sp_max = sp
  next_Face_id = start_Face

  while searching
    searching = false # abort condition
    Face = Grid.Faces[start_Face]
    for i = 1 : length(Face.Stencil)
      Face_id = Face.Stencil[i]  
      if Face_id != start_Face                                   # 0 is the "undefined" value for the face id
        neighbour_cc = Grid.Faces[Face_id].Mid       # get cartesian coordinates of neighbour face
        sp = dot(target_cc, neighbour_cc)         # calculate measure for distance to target point
        if sp > sp_max
          sp_max = sp                                # save new distance measure
          next_Face_id = Face_id                     # save face id
          searching = true                           # continue with search loop
        end  
      end  
    end  
    start_Face = next_Face_id
  end
  Inside = InsideFace(target_cc,Grid.Faces[start_Face],Grid)
  if Inside == false
    for i = 1 : length(Grid.Faces[start_Face].Stencil)
      iF = Grid.Faces[start_Face].Stencil[i]  
      Inside = InsideFace(target_cc,Grid.Faces[iF],Grid)
      if Inside == true
        start_Face = iF  
        break  
      end  
    end  
  end  
  (iPosFace, jPosFace) = walk_to_face(target_cc,Grid.Faces[start_Face],Trans,xw,Rad)
  return (start_Face, iPosFace, jPosFace)       # move one face toward target point
end 

function walk_to_face(target_cc,Face,Trans,xw,Rad)

  iPos = 1
  jPos = 1
  P = Point(Trans(xw[1],xw[1],Face,Rad))
  sp_max = dot(target_cc, P)
  for j = 1 : length(xw)
    for i = 1 : length(xw)
      P = Point(Trans(xw[i],xw[j],Face,Rad))
      sp = dot(target_cc, P)
      if sp > sp_max
        iPos = i
        jPos = j
        sp_max = sp
      end
    end
  end
  return (iPos, jPos)
end
