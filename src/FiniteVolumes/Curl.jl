function CurlNode2!(c,uN,Metric,Grid)

  FT = eltype(c)
  t = zeros(FT,3)
  for iE = 1 : Grid.NumEdges
    N1 =  Grid.Edges[iE].N[1]
    c[N1] += Metric.DualEdge[iE]*uN[iE]
    N2 =  Grid.Edges[iE].N[2]
    c[N2] -= Metric.DualEdge[iE]*uN[iE]
  end
  @. c /= Metric.DualVolume
end

function CurlNode1!(c,uN,Metric,Grid)

# u = utL * tL + utR * tR
# (u,nL) = utR * (tR,nL)
# utR = unL / (tR,nL)
# (u,nR) = utL * (tL,nR)
# utL = unR / (tL,nR)
# (u,s) = utL * (tL,s) + utR * (tR,s)
#       = unR / (tL,nR) * (tL,s) + unL / (tR,nL) * (tR,s)

  Rad = Grid.Rad
  FT = eltype(c)
  t = zeros(FT,3)
  @. c = 0
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      iER= Grid.Faces[iF].E[i]
      if i == 1
        iEL = Grid.Faces[iF].E[end]
      else
        iEL = Grid.Faces[iF].E[i-1]
      end

      uuL = uN[iEL]
      nL = Grid.Edges[iEL].n
      tL = Grid.Edges[iEL].t
      sL = Grid.Faces[iF].Mid - Grid.Edges[iEL].Mid
      sL = sL / norm(sL) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iEL].Mid) * Grid.Rad
      uuR = uN[iER]
      nR = Grid.Edges[iER].n
      tR = Grid.Edges[iER].t
      sR = Grid.Faces[iF].Mid - Grid.Edges[iER].Mid
      sR = sR / norm(sR) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iER].Mid) * Grid.Rad

      usR = uuR / dot(nR,tL) * dot(tL,sR) + uuL / dot(nL,tR) * dot(tR,sR)  
      usL = uuR / dot(nR,tL) * dot(tL,sL) + uuL / dot(nL,tR) * dot(tR,sL)  

      c[iN] += usR - usL
    end
  end
  @. c /= Metric.DualVolume
end
function CurlNode!(c,uN,Metric,Grid)

# u = utL * tL + utR * tR
# (u,nL) = utR * (tR,nL)
# utR = unL / (tR,nL)
# (u,nR) = utL * (tL,nR)
# utL = unR / (tL,nR)
# (u,s) = utL * (tL,s) + utR * (tR,s)
#       = unR / (tL,nR) * (tL,s) + unL / (tR,nL) * (tR,s)

  Rad = Grid.Rad
  FT = eltype(c)
  t = zeros(FT,3)
  @. c = 0
  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      iER= Grid.Faces[iF].E[i]
      if i == 1
        iEL = Grid.Faces[iF].E[end]
      else
        iEL = Grid.Faces[iF].E[i-1]
      end
      if i == 1
        iNL = Grid.Faces[iF].N[end]
        iNR = Grid.Faces[iF].N[i+1]
      elseif i == length(Grid.Faces[iF].N)
        iNL = Grid.Faces[iF].N[i-1]
        iNR = Grid.Faces[iF].N[1]
      else
        iNL = Grid.Faces[iF].N[i-1]
        iNR = Grid.Faces[iF].N[i+1]
      end    

      uuL = uN[iEL]
      nL = Grid.Edges[iEL].n
      tL = Grid.Edges[iEL].t
      sL = Grid.Faces[iF].Mid - Grid.Edges[iEL].Mid
      sL = sL / norm(sL) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iEL].Mid) * Grid.Rad
      uuR = uN[iER]
      nR = Grid.Edges[iER].n
      tR = Grid.Edges[iER].t
      sR = Grid.Faces[iF].Mid - Grid.Edges[iER].Mid
      sR = sR / norm(sR) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iER].Mid) * Grid.Rad

      usR = uuR / dot(nR,tL) * dot(tL,sR) + uuL / dot(nL,tR) * dot(tR,sR)  
      usL = uuR / dot(nR,tL) * dot(tL,sL) + uuL / dot(nL,tR) * dot(tR,sL)  

      c[iN] += 0.5 * (usR - usL)
      c[iNL] += 0.5 * usL
      c[iNR] += -0.5 * usR 

    end
  end
  @. c /= Metric.DualVolume
end
function CurlNodeMatrix(Metric,Grid)

# u = utL * tL + utR * tR
# (u,nL) = utR * (tR,nL)
# utR = unL / (tR,nL)
# (u,nR) = utL * (tL,nR)
# utL = unR / (tL,nR)
# (u,s) = utL * (tL,s) + utR * (tR,s)
#       = unR / (tL,nR) * (tL,s) + unL / (tR,nL) * (tR,s)

  Rad = Grid.Rad
  DV = Metric.DualVolume
  FT = eltype(DV)
  t = zeros(FT,3)
  RowInd = Int64[]
  ColInd = Int64[]
  Val = FT[]

  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      iER= Grid.Faces[iF].E[i]
      if i == 1
        iEL = Grid.Faces[iF].E[end]
      else
        iEL = Grid.Faces[iF].E[i-1]
      end
      if i == 1
        iNL = Grid.Faces[iF].N[end]
        iNR = Grid.Faces[iF].N[i+1]
      elseif i == length(Grid.Faces[iF].N)
        iNL = Grid.Faces[iF].N[i-1]
        iNR = Grid.Faces[iF].N[1]
      else
        iNL = Grid.Faces[iF].N[i-1]
        iNR = Grid.Faces[iF].N[i+1]
      end    

      nL = Grid.Edges[iEL].n
      tL = Grid.Edges[iEL].t
      sL = Grid.Faces[iF].Mid - Grid.Edges[iEL].Mid
      sL = sL / norm(sL) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iEL].Mid) * Grid.Rad
      nR = Grid.Edges[iER].n
      tR = Grid.Edges[iER].t
      sR = Grid.Faces[iF].Mid - Grid.Edges[iER].Mid
      sR = sR / norm(sR) * Grids.SizeGreatCircle(Grid.Faces[iF].Mid,Grid.Edges[iER].Mid) * Grid.Rad

      pRR = dot(tL,sR) / dot(nR,tL)
      pRL = dot(tR,sR) / dot(nL,tR)
      pLR = dot(tL,sL) / dot(nR,tL)
      pLL = dot(tR,sL) / dot(nL,tR)

#     usR = uuR / dot(nR,tL) * dot(tL,sR) + uuL / dot(nL,tR) * dot(tR,sR)  
#     usL = uuR / dot(nR,tL) * dot(tL,sL) + uuL / dot(nL,tR) * dot(tR,sL)  

#     c[iN] += 0.5 * (usR - usL)
#     c[iNL] += 0.5 * usL
#     c[iNR] += -0.5 * usR 
      push!(RowInd,iN)
      push!(ColInd,iER)
      push!(Val,0.5*pRR/DV[iN])
      push!(RowInd,iN)
      push!(ColInd,iEL)
      push!(Val,0.5*pRL/DV[iN])

      push!(RowInd,iN)
      push!(ColInd,iER)
      push!(Val,-0.5*pLR/DV[iN])
      push!(RowInd,iN)
      push!(ColInd,iEL)
      push!(Val,-0.5*pLL/DV[iN])

      push!(RowInd,iNR)
      push!(ColInd,iER)
      push!(Val,-0.5*pRR/DV[iNR])
      push!(RowInd,iNR)
      push!(ColInd,iEL)
      push!(Val,-0.5*pRL/DV[iNR])

      push!(RowInd,iNL)
      push!(ColInd,iER)
      push!(Val,0.5*pLR/DV[iNL])
      push!(RowInd,iNL)
      push!(ColInd,iEL)
      push!(Val,0.5*pLL/DV[iNL])
    end
  end

  M = sparse(RowInd, ColInd, Val)
end

function TagentialVelocity(uT,uN,Metric,Grid)

# u = unL * nL + unR * nR
# (u,tL) = utR * (tL,nR)
# (u,tR) = utL * (tR,nL)

# (u,tL) = utL * (tL,tL) + utR * (tR,tL)
# (u,tR) = utL * (tL,tR) + utR * (tR,tR)

  FT = eltype(Metric.DualVolume)
  @. uT = 0

  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      iER= Grid.Faces[iF].E[i]
      if i == 1
        iEL = Grid.Faces[iF].E[end]
      else
        iEL = Grid.Faces[iF].E[i-1]
      end

      nL = Grid.Edges[iEL].n
      tL = Grid.Edges[iEL].t
      uuL = uN[iEL]
      nR = Grid.Edges[iER].n
      tR = Grid.Edges[iER].t
      uuR = uN[iER]

      uTR = uuR / dot(nR,tL) * dot(tL,tR) + uuL / dot(nL,tR) * dot(tR,tR)  
      uTL = uuR / dot(nR,tL) * dot(tL,tL) + uuL / dot(nL,tR) * dot(tR,tL)  

      uT[iER] += 0.25 * uTR
      uT[iEL] += 0.25 * uTL
    end
  end
end

function TagentialVelocity2(uT,uN,Metric,Grid)


  FT = eltype(Metric.DualVolume)
  Rad = Grid.Rad
  @. uT = 0

  for iF = 1 : Grid.NumFaces
    for iE in Grid.Faces[iF].E
      u = Grids.Point(0.0,0.0,0.0)  
      Mid = Grid.Edges[iE].Mid
      for j = 1 : length(Grid.Faces[iF].E)
        jE = Grid.Faces[iF].E[j]  
        if jE != iE  
          t = (Grid.Edges[jE].Mid - Mid)  
          t = t / norm(t) * Grids.SizeGreatCircle(Mid,Grid.Edges[jE].Mid) * Grid.Rad
          u = u - (uN[jE] * Metric.PrimalEdge[jE] * Grid.Faces[iF].OrientE[j]) * t
        end
      end
      u /= Metric.PrimalVolume[iF]
      k = Mid / norm(Mid)
      u = u - dot(u,k) * k
      uTE = dot(Grid.Edges[iE].t,u)
      uT[iE] += 0.5 * uTE 
    end  
  end    
end

function TagentialVelocity1(uT,uN,Metric,Grid)

# u = unL * nL + unR * nR
# (u,tL) = utR * (tL,nR)
# (u,tR) = utL * (tR,nL)

# (u,tL) = utL * (tL,tL) + utR * (tR,tL)
# (u,tR) = utL * (tL,tR) + utR * (tR,tR)

  FT = eltype(Metric.DualVolume)
  Rad = Grid.Rad
  @. uT = 0

  for iF = 1 : Grid.NumFaces
    MidF = Grid.Faces[iF].Mid  
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      P = Grid.Nodes[iN].P
      iER= Grid.Faces[iF].E[i]
      if i == 1
        iEL = Grid.Faces[iF].E[end]
      else
        iEL = Grid.Faces[iF].E[i-1]
      end

      k = Grid.Nodes[iN].P
      iN1 = Grid.Edges[iEL].N[1]
      iN2 = Grid.Edges[iEL].N[2]
      t = Grid.Nodes[iN2].P - Grid.Nodes[iN1].P
      n = Grids.cross(k,t)
      t = Grids.cross(n,k)
      nL = n / Grids.norm(n)
      tL = t / Grids.norm(t)
      MidEL = Grid.Edges[iEL].Mid
      EdgeELVol = Grids.AreaSphericalTriangle(MidEL,MidF,P) * Rad^2
      iN1 = Grid.Edges[iER].N[1]
      iN2 = Grid.Edges[iER].N[2]
      t = Grid.Nodes[iN2].P - Grid.Nodes[iN1].P
      n = Grids.cross(k,t)
      t = Grids.cross(n,k)
      nR = n / Grids.norm(n)
      tR = t / Grids.norm(t)
      MidER = Grid.Edges[iER].Mid
      EdgeERVol = Grids.AreaSphericalTriangle(MidER,MidF,P) * Rad^2

      uuL = uN[iEL]
      uuR = uN[iER]

      uTR = uuL / dot(nL,tR) * dot(tR,Grid.Edges[iER].t) 
      uTL = uuR / dot(nR,tL) * dot(tL,Grid.Edges[iEL].t)

      uT[iER] += EdgeERVol * uTR 
      uT[iEL] += EdgeELVol * uTL
    end
  end
  @. @views uT /= (Metric.DualEdgeVolume[1,:] + Metric.DualEdgeVolume[2,:])
end

function TagentialVelocityMatrix(Metric,Grid)

# u = utL * tL + utR * tR
# (u,nL) = utR * (tR,nL)
# utR = unL / (tR,nL)
# (u,nR) = utL * (tL,nR)
# utL = unR / (tL,nR)

  FT = eltype(Metric.DualVolume)
  t = zeros(FT,3)
  RowInd = Int64[]
  ColInd = Int64[]
  Val = FT[]

  for iF = 1 : Grid.NumFaces
    for i = 1 : length(Grid.Faces[iF].N)
      iN = Grid.Faces[iF].N[i]  
      iER= Grid.Faces[iF].E[i]
      if i == 1
        iEL = Grid.Faces[iF].E[end]
      else
        iEL = Grid.Faces[iF].E[i-1]
      end

      nL = Grid.Edges[iEL].n
      tL = Grid.Edges[iEL].t
      nR = Grid.Edges[iER].n
      tR = Grid.Edges[iER].t

      utL = 1.0 / dot(nR,tL)
      utR = 1.0 / dot(nL,tR)

#     uTR = uuR / dot(nR,tL) * dot(tL,tR) + uuL / dot(nL,tR) * dot(tR,tR)  
#     uTL = uuR / dot(nR,tL) * dot(tL,tL) + uuL / dot(nL,tR) * dot(tR,tL)  

#     uI[iER] += 0.25 * utR
#     uI[iEL] += 0.25 * utL
      push!(RowInd,iER)
      push!(ColInd,iEL)
      temp = 0.25 / dot(nL,tR) 
      push!(Val,temp)
      push!(RowInd,iER)
      push!(ColInd,iER)
      temp = 0.25 * dot(tL,tR) / dot(nR,tL) 
      push!(Val,temp)

      push!(RowInd,iEL)
      push!(ColInd,iEL)
      temp = 0.25 * dot(tR,tL) / dot(nL,tR)
      push!(Val,temp)
      push!(RowInd,iEL)
      push!(ColInd,iER)
      temp = 0.25 / dot(nR,tL)
      push!(Val,temp)
    end
  end

  sparse(RowInd, ColInd, Val)
end

function TagentialVelocity2Matrix(Metric,Grid)


  FT = eltype(Metric.DualVolume)
  Rad = Grid.Rad

  RowInd = Int64[]
  ColInd = Int64[]
  Val = FT[]
  for iF = 1 : Grid.NumFaces
    uN = zeros(length(Grid.Faces[iF].E))  
    for iE in Grid.Faces[iF].E
      u = Grids.Point(0.0,0.0,0.0)
      Mid = Grid.Edges[iE].Mid
      for k = 1 : length(Grid.Faces[iF].E)
        kE = Grid.Faces[iF].E[k]  
        @. uN = 0  
        uN[k] = 1
        for j = 1 : length(Grid.Faces[iF].E)
          jE = Grid.Faces[iF].E[j]
          if jE != iE
            t = (Grid.Edges[jE].Mid - Mid)
            t = t / norm(t) * Grids.SizeGreatCircle(Mid,Grid.Edges[jE].Mid) * Grid.Rad
            u = u - (uN[j] * Metric.PrimalEdge[jE] * Grid.Faces[iF].OrientE[j]) * t
          end
        end
        u /= Metric.PrimalVolume[iF]
        k = Mid / norm(Mid)
        u = u - dot(u,k) * k
        uTE = 0.5 * dot(Grid.Edges[iE].t,u)
        push!(RowInd,iE)
        push!(ColInd,kE)
        push!(Val,uTE)
      end  
    end
  end
  sparse(RowInd, ColInd, Val)
end
