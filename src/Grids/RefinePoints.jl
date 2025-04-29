function PrintPoints(Refine,::Grids.Quad)

  NumRefine = (Refine + 1) * (Refine + 1)
  RefineMidPoints = zeros(NumRefine,2)
  RefinePoints = zeros(NumRefine,4,2)
  dksi = 2 / (Refine + 1)
  iRefine = 0
  for jR in 1 : Refine + 1
    for iR in 1 : Refine + 1
      iRefine += 1
      RefinePoints[iRefine,1,1] = -1.0 + (iR - 1) * dksi
      RefinePoints[iRefine,1,2] = -1.0 + (jR - 1) * dksi
      RefinePoints[iRefine,2,1] = -1.0 + (iR    ) * dksi
      RefinePoints[iRefine,2,2] = -1.0 + (jR - 1) * dksi
      RefinePoints[iRefine,3,1] = -1.0 + (iR    ) * dksi
      RefinePoints[iRefine,3,2] = -1.0 + (jR    ) * dksi
      RefinePoints[iRefine,4,1] = -1.0 + (iR - 1) * dksi
      RefinePoints[iRefine,4,2] = -1.0 + (jR    ) * dksi
      RefineMidPoints[iRefine,1] = 0.25 * (RefinePoints[iRefine,1,1] +
        RefinePoints[iRefine,2,1] + RefinePoints[iRefine,3,1] +
        RefinePoints[iRefine,4,1])
      RefineMidPoints[iRefine,2] = 0.25 * (RefinePoints[iRefine,1,2] +
        RefinePoints[iRefine,2,2] + RefinePoints[iRefine,3,2] +
        RefinePoints[iRefine,4,2])
    end
  end
  return RefinePoints, RefineMidPoints
end  

function PrintPoints(Refine,::Grids.Tri)

  NumRefine = (Refine + 1) * (Refine + 1)
  RefineMidPoints = zeros(NumRefine,2)
  RefinePoints = zeros(NumRefine,3,2)
  dksi = 2 / (Refine + 1)
  iRefine = 0
  for jRR in Refine + 1 : -1 : 1
    jR = Refine + 2 - jRR
    for iR in 1 : jRR
      iRefine += 1
      RefinePoints[iRefine,1,1] = -1.0 + (iR - 1) * dksi
      RefinePoints[iRefine,1,2] = -1.0 + (jR - 1) * dksi
      RefinePoints[iRefine,2,1] = -1.0 + (iR    ) * dksi
      RefinePoints[iRefine,2,2] = -1.0 + (jR - 1) * dksi
      RefinePoints[iRefine,3,1] = -1.0 + (iR - 1) * dksi
      RefinePoints[iRefine,3,2] = -1.0 + (jR    ) * dksi
      RefineMidPoints[iRefine,1] = 1/3 * (RefinePoints[iRefine,1,1] +
        RefinePoints[iRefine,2,1] + RefinePoints[iRefine,3,1])
      RefineMidPoints[iRefine,2] = 1/3 * (RefinePoints[iRefine,1,2] +
        RefinePoints[iRefine,2,2] + RefinePoints[iRefine,3,2])
      if iR == jRR
        continue
      end
      iRefine += 1
      RefinePoints[iRefine,1,1] = -1.0 + (iR    ) * dksi
      RefinePoints[iRefine,1,2] = -1.0 + (jR - 1) * dksi
      RefinePoints[iRefine,2,1] = -1.0 + (iR    ) * dksi
      RefinePoints[iRefine,2,2] = -1.0 + (jR    ) * dksi
      RefinePoints[iRefine,3,1] = -1.0 + (iR - 1) * dksi
      RefinePoints[iRefine,3,2] = -1.0 + (jR    ) * dksi
      RefineMidPoints[iRefine,1] = 1/3 * (RefinePoints[iRefine,1,1] +
        RefinePoints[iRefine,2,1] + RefinePoints[iRefine,3,1])
      RefineMidPoints[iRefine,2] = 1/3 * (RefinePoints[iRefine,1,2] +
        RefinePoints[iRefine,2,2] + RefinePoints[iRefine,3,2])
    end
  end
  return RefinePoints, RefineMidPoints
end

