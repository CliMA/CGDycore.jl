function MetricFace(JdXdxI,NormalDir,TangDir)

#
#  JdXdxI(1:OrdPoly+1,1:OrdPoly+1,3,3)
#

  for j = 1 : OrdPoly + 1
    for i = 1 : OrdPoly + 1
      SurfElem[i,j] = norm(JdXdxI[i,j,NormalDir,:])
      NormVec[i,j,:] = JdXdxI[i,j,NormalDir,:] / SurfElem[i,j]
      Tang1Vec[i,j,:] = JdXdxI[i,j,TangDir,:] - dot(JdXdxI[i,j,TangDir,:],NormVec[i,j,:]) * NormVec[i,j,:] 
      Tang1Vec[i,j,:] = Tang1Vec[i,j,:] / norm(Tang1Vec[i,j,:])
      Tang2Vec[i,j,:] = cross(NormVec[i,j,:], Tang1Vec[i,j,:])
    end   
  end   
  return (SurfElem, NormVec, Tang1Vec, Tang2Vec)
end   
