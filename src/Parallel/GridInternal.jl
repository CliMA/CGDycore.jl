  function CreateCubedSphereMesh_Internal!(dm::DM, R::Float64)
    numCells = 6
    numEdges = 12
    numVerts = 8
    firstVertex = numCells
    firstEdge = numCells + numVerts
    # Build Topology 
    SetChart!(dm, 0, numCells+numEdges+numVerts)
    for c = 1 : numCells
      SetConeSize!(dm, c, 4)
    end  
    for e = firstEdge + 1 : firstEdge+numEdges
      SetConeSize!(dm, e, 2)
    end  
    #DMSetUp(dm)
    DMSetUp_Plex!(dm)
    cone=zeros(Int,4)
    ornt=zeros(Int,4)

    # Cell 1 */
    cone[1] = 15
    cone[2] = 16
    cone[3] = 17
    cone[4] = 18
    SetCone!(dm, 1, cone)
    ornt[1] = 0
    ornt[2] = 0
    ornt[3] = 0
    ornt[4] = 0
    SetConeOrientation!(dm, 1, ornt)

    # Cell 2 */
    cone[1] = 19
    cone[2] = 20
    cone[3] = 15
    cone[4] = 21
    SetCone!(dm, 2, cone)
    ornt[1] = 0
    ornt[2] = 0
    ornt[3] = -1
    ornt[4] = 0
    SetConeOrientation!(dm, 2, ornt)

    #Cell 3 */
    cone[1] = 22
    cone[2] = 23
    cone[3] = 19
    cone[4] = 24
    SetCone!(dm, 3, cone)
    ornt[1] = 0
    ornt[2] = 0
    ornt[3] = -1
    ornt[4] = 0
    SetConeOrientation!(dm, 3, ornt)

    #Cell 4 */
    cone[1] = 20
    cone[2] = 23
    cone[3] = 25
    cone[4] = 16
    SetCone!(dm, 4, cone)
    ornt[1] = -1
    ornt[2] = -1
    ornt[3] = 0
    ornt[4] = -1
    SetConeOrientation!(dm, 4, ornt)

    #Cell 5 */
    cone[1] = 17
    cone[2] = 25
    cone[3] = 22
    cone[4] = 26
    SetCone!(dm, 5, cone)
    ornt[1] = -1
    ornt[1] = -1
    ornt[3] = -1
    ornt[4] = 0
    SetConeOrientation!(dm, 5, ornt)

    #Cell 6 */
    cone[1] = 21
    cone[2] = 18
    cone[3] = 26
    cone[4] = 24
    SetCone!(dm, 6, cone)
    ornt[1] = -1
    ornt[2] = -1
    ornt[3] = -1
    ornt[4] = -1
    SetConeOrientation!(dm, 6, ornt)

    # Edges */
    cone[1] =  7; cone[2] =  8;
    SetCone!(dm, 15, cone)
    cone[1] =  8; cone[2] =  9;
    SetCone!(dm, 16, cone)
    cone[1] =  9; cone[2] =  10;
    SetCone!(dm, 17, cone)
    cone[1] =  10; cone[2] =  7;
    SetCone!(dm, 18, cone)
    cone[1] = 11; cone[2] = 12;
    SetCone!(dm, 19, cone)
    cone[1] = 12; cone[2] =  8;
    SetCone!(dm, 20, cone)
    cone[1] =  7; cone[2] = 11;
    SetCone!(dm, 21, cone)
    cone[1] = 13; cone[2] = 14;
    SetCone!(dm, 22, cone)
    cone[1] = 14; cone[2] = 12;
    SetCone!(dm, 23, cone)
    cone[1] = 11; cone[2] = 13;
    SetCone!(dm, 24, cone)
    cone[1] = 14; cone[2] =  9;
    SetCone!(dm, 25, cone)
    cone[1] = 13; cone[2] =  10;
    SetCone!(dm, 26, cone)
    Symmetrize(dm)
    Stratify(dm)
  end
