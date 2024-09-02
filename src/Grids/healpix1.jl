using LinearAlgebra

function HealpixGrid(ns)
  """
  healpix.jl: Generator of HEALPix grid

  Healpix points are plotted on the sphere or in the lon-lat plane.
  Rhomboids are created with their centroids in healpix points.
  On the sphere, the edges of the rhomboids are great arcs.
  """

  # Parameters
  #ns = 2  # for healpix H<ns>

  pi = π

  npix = 12 * ns * ns
  nrings = 4 * ns - 1
  rad2deg = 180 / pi
  deg2rad = 1 / rad2deg
  ncap_ids = zeros(Int, 2)

  # HEALPix dual points
  dualp = zeros(npix, 2)

  # North cap
  for p in 1:npix
    ph = 0.5 * p
    i = floor(Int,sqrt(ph - sqrt(floor(ph)))) + 1
    if i > ns
      ncap_ids[1] = p - 1  # Adjust for Julia's 1-based index
      break
    end
    j = p - 2 * i * (i - 1)
    dualp[p, 1] = acos(1 - i^2 / (3 * ns^2)) * rad2deg
    dualp[p, 2] = (0.25 * pi * (2 * j - 1) / i) * rad2deg
    @show "10",p,dualp[p, :]
  end

  # North belt
  for p in ncap_ids[1] + 1:npix
    ph = p - 1 - 2 * ns * (ns - 1)
    i = floor(Int,0.25 * ph / ns) + ns
    if i > 2 * ns
      ncap_ids[2] = p - 1  # Adjust for Julia's 1-based index
      break
    end
    j = 1 + ph % (4 * ns)
    s = (i - ns + 1) % 2
    dualp[p, 1] = acos((4 * ns - 2 * i) / (3 * ns)) * rad2deg
    dualp[p, 2] = (0.25 * pi * (2 * j - s) / ns) * rad2deg
    @show "20",p,dualp[p, :],ph,i,j,s
  end

  # South belt & south cap
  for p in ncap_ids[2] + 1:npix
    dualp[p, 1] = 180 - dualp[npix - p + 1, 1]
    dualp[p, 2] = dualp[npix - p + 1, 2]
  end
  for p = 1 : size(dualp,1)
    @show p,dualp[p,:]  
  end
  stop

  # HEALPix vertex points
  vertp = zeros(npix + 2, 2)
  vertpC = zeros(npix + 2, 3)

  vertp[1, 1] = 0
  vertp[1, 2] = 180
  vertp[npix + 2, 1] = 180
  vertp[npix + 2, 2] = 180

  # North cap
  for p in 1:ncap_ids[1]
    ph = 0.5 * p
    i = floor(Int,sqrt(ph - sqrt(floor(ph)))) + 1
    last_id = 2 * i * (i - 1)
    j = p - last_id
    vertp[p + 1, 1] = dualp[p, 1]
    pp = p > 1 ? p - 1 : last_id + 4 * i - 1
    vertp[p + 1, 2] = dualp[pp, 2] + dualp[p, 2]
    vertp[p + 1, 2] = p > pp ? 0.5 * vertp[p + 1, 2] : 0.5 * (vertp[p + 1, 2] + 360)
    @show p+1,vertp[p+1,:]
  end

  # North belt
  for p in ncap_ids[1] + 1:ncap_ids[2]
    ph = p - 2 * ns * (ns - 1)
    i = floor(Int,0.25 * ph / ns) + ns
    j = 1 + ph % (4 * ns)
    vertp[p + 1, 1] = dualp[p, 1]
    pp = p > 1 ? p - 1 : ncap_ids[1] + 4 * ns * (i - ns) - 1
    vertp[p + 1, 2] = dualp[pp, 2] + dualp[p, 2]
    vertp[p + 1, 2] = p > pp ? 0.5 * vertp[p + 1, 2] : 0.5 * (vertp[p + 1, 2] + 360)
    @show p+1,vertp[p+1,:],dualp[p,1]
  end

  # South belt & south cap
  for p in ncap_ids[2] + 1:npix
    vertp[p + 1, 1] = 180 - vertp[npix - p + 1, 1]
    vertp[p + 1, 2] = vertp[npix - p + 1, 2]
  end

  for p = 1 : size(vertp,1)
    lon = vertp[p,1]
    lat = vertp[p,2]
    vertpC[p,1] = sind(lon)*cosd(lat)
    vertpC[p,2] = cosd(lon)*cosd(lat)
    vertpC[p,3] = sind(lat)
    @show p,vertp[p,:]
  end    
  stop

  # HEALPix cells
  el_vid = zeros(Int, npix, 4)

  # North cap for cells
  for el in 1:ncap_ids[1]
    ph = 0.5 * el
    i = floor(sqrt(ph - sqrt(floor(ph)))) + 1
    j = el - 2 * i * (i - 1)
    if i > 1
      if j < 4 * i
        el_vid[el,1] = j-floor(Int,((j-1)/i))
      else
        el_vid[el,1] = 1
      end
    else
      el_vid[el,1] = 0
    end
    if i > 2
      el_vid[el,1] += 2*(i-1)*(i-2)
    end
    el_vid[el, 2] = j + 2 * i * (i - 1)
    el_vid[el, 4] = j < 4 * i ? j + 1 + 2 * i * (i - 1) : 1 + 2 * i * (i - 1)
#   el_vid[el, 3] = el < 2 * (ns - 1) * ns ? 2 * (i + 1) * i + j + 1 : 2 * (i + 1) * i + j
    if el - 1 < 2*(ns-1)*ns 
        el_vid[el,3] = 2*(i+1)*i + j + 1 + floor(Int,(j-1)/i)
    else
        el_vid[el,3] = 2*(i+1)*i + j
    end    
  end

  # North belt for cells
  for el in ncap_ids[1] + 1:6 * ns * ns + 2 * ns
    ph = el - 1 - 2 * ns * (ns - 1)
    i = floor(0.25 * ph / ns) + ns
    j = 1 + ph % (4 * ns)
    s = (i - ns + 1) % 2
    el_vid[el, 1] = 2 * ns * (ns - 1) + 4 * ns * (i - ns - 1) + 1
    if j == 4 * ns && s == 0
      el_vid[el, 1] -= s  
    else
      el_vid[el, 1] += j - s  
    end  
    el_vid[el, 2] = el - 1 + 1
    el_vid[el, 4] = el - 1 + 2 - (j == 4 * ns ? 4 * ns : 0)
    if el - 1 < 6*ns*ns-2*ns 
      el_vid[el,3] = 1 - s 
      if j == 4 * ns && s == 0
        el_vid[el,3] += 2*ns*(ns-1)+4*ns*(i-ns+1) 
      else
         el_vid[el,3] += el - 1+4*ns+1 
      end   
    end  
  end

# Equator for cells
  for el in 6 * ns * ns - 2 * ns + 1:6 * ns * ns + 2 * ns
    el_vid[el, 3] = npix + 1 - el_vid[el, 1]
  end

# South belt & south cap for cells
  for el in 6 * ns * ns + 2 * ns + 1:npix
    el_vid[el, 1] = npix + 1 - el_vid[npix - el + 1, 1]
    el_vid[el, 2] = npix + 1 - el_vid[npix - el + 1, 2]
    el_vid[el, 4] = npix + 1 - el_vid[npix - el + 1, 4]
#   el_vid[el, 3] = el_vid[npix - el + 1, 3] > 6 * ns * ns - 2 * ns && el_vid[npix - el + 1, 3] <= 
#     6 * ns * ns + 2 * ns ? el_vid[npix - el + 1, 3] : npix + 1 - el_vid[npix - el + 1, 3]
    if el_vid[npix-el+1,3] > 6*ns*ns-2*ns && el_vid[npix-el+1,3] <= 6*ns*ns+2*ns 
        el_vid[el,3] = el_vid[npix-el+1,3]
    else
        el_vid[el,3] = npix + 1 - el_vid[npix-el+1,3]  
    end    
  end
# for el = 1 : size(el_vid,1)
#    @show el,el_vid[el,:]
# end 

  @. el_vid += 1
  stop

  Node2Edge = Dict()
  E1 = zeros(Int,0)
  E2 = zeros(Int,0)
  iE = 1
  for el = 1 : size(el_vid,1)
    # Edge 1
    n1 = el_vid[el,1]
    n2 = el_vid[el,2]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
      push!(E1,n1)
      push!(E2,n2)
    end  
    # Edge 2
    n1 = el_vid[el,2]
    n2 = el_vid[el,3]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
    end  
    # Edge 3
    n1 = el_vid[el,3]
    n2 = el_vid[el,4]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
    end  
    # Edge 4
    n1 = el_vid[el,4]
    n2 = el_vid[el,1]
    if n1 < n2
      Node2Edge[[n1,n2]] = iE
      iE += 1
    end  
  end   
  F = zeros(Int,size(el_vid))
  for el = 1 : size(el_vid,1)
    # Edge 1
    n1 = el_vid[el,1]
    n2 = el_vid[el,2]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    @show "V1",el
    F[el,1] = Node2Edge[e]
    # Edge 2
    n1 = el_vid[el,2]
    n2 = el_vid[el,3]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    @show "V2",el
    F[el,2] = Node2Edge[e]
    # Edge 3
    n1 = el_vid[el,3]
    n2 = el_vid[el,4]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    @show "V3",el
    F[el,3] = Node2Edge[e]
    # Edge 4
    n1 = el_vid[el,4]
    n2 = el_vid[el,1]
    if n1 < n2
      e = [n1,n2]  
    else
      e = [n2,n1]  
    end  
    @show "V4",el
    F[el,4] = Node2Edge[e]
  end

end
