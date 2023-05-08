function area_of_cap(s_cap::Float64) 
    #
    # AREA_OF_CAP Area of spherical cap
    #
    # AREA_OF_CAP(S_CAP) sets AREA to be the area of an S^2 spherical
    # cap of spherical radius S_CAP.
    #
    return 4.0 * pi * sin(0.5 * s_cap)^2
end    

function area_of_collar(a_top::Float64,a_bot::Float64) 
    # AREA_OF_COLLAR Area of spherical collar
    #
    # AREA_OF_COLLAR(A_TOP, A_BOT) sets AREA to be the area of an S^2 spherical
    # collar specified by A_TOP, A_BOT, where A_TOP is top (smaller) spherical
    # radius,
    # A_BOT is bottom (larger) spherical radius.
    #
    return area_of_cap(a_bot) - area_of_cap(a_top)
end

function sradius_of_cap(area::Float64) 
    # SRADIUS_OF_CAP(AREA) returns the spherical radius of
    # an S^2 spherical cap of area AREA.
    #
    return 2. * asin(0.5 * sqrt(area / pi))
end

function area_of_ideal_region(N::Int) 
    #
    # AREA_OF_IDEAL_REGION(N) sets AREA to be the area of one of N equal
    # area regions on S^2, that is 1/N times AREA_OF_SPHERE.
    #
    area_of_sphere = 4*pi
    return area_of_sphere / N
end

function polar_colat(N::Int) 
    #
    # Given N, determine the colatitude of the North polar spherical cap.
    #
    if N == 1
        polar_c = pi
    end
    if N == 2
        polar_c = 0.5 * pi
    end
    if N > 2
        polar_c = sradius_of_cap(area_of_ideal_region(N))
    end
    return polar_c
end

function ideal_collar_angle(N::Int) 
    #
    # IDEAL_COLLAR_ANGLE The ideal angle for spherical collars of an EQ partition
    #
    # IDEAL_COLLAR_ANGLE(N) sets ANGLE to the ideal angle for the
    # spherical collars of an EQ partition of the unit sphere S^2 into N regions.
    #
    return sqrt(area_of_ideal_region(N))
end

function ideal_region_list!(N::Int,c_polar::Float64,n_collars::Int,r_regions::Array{Float64,1})
    #
    # IDEAL_REGION_LIST The ideal real number of regions in each zone
    #
    #  List the ideal real number of regions in each collar, plus the polar caps.
    #
    #  Given N, c_polar and n_collars, determine r_regions, a list of the ideal
    #  real
    #  number of regions in each collar, plus the polar caps.
    #  The number of elements is n_collars+2.
    #  r_regions[1] is 1.
    #  r_regions[n_collars+2] is 1.
    #  The sum of r_regions is N.
    #
    # real(wp),intent(out) :: r_regions(n_collars+2)
    r_regions[1] = 1.
    if n_collars > 0
        #
        # Based on n_collars and c_polar, determine a_fitting,
        # the collar angle such that n_collars collars fit between the polar caps.
        #
        a_fitting  = (pi - 2. * c_polar) / n_collars
        ideal_region_area = area_of_ideal_region(N)
        for collar_n = 1:n_collars
            ideal_collar_area = area_of_collar(c_polar + (collar_n - 1) * a_fitting, c_polar + collar_n * a_fitting)
            r_regions[collar_n+1] = ideal_collar_area / ideal_region_area
        end
    end
    r_regions[2 + n_collars] = 1.
end

function num_collars(N::Int, c_polar::Float64, a_ideal::Float64)
    #
    # NUM_COLLARS The number of collars between the polar caps
    #
    #  Given N, an ideal angle, and c_polar,
    #  determine n_collars, the number of collars between the polar caps.
    #
    enough = (N > 2) && (a_ideal > 0)
    if enough
        return max(1, round(Int,(pi - 2. * c_polar) / a_ideal))
    else  
        return 1
    end
end

function round_to_naturals!(N::Int, ncollars::Int, r_regions::Array{Float64,1}, n_regions::Array{Int,1})
    # ROUND_TO_NATURALS Round off a given list of numbers of regions
    #
    #  Given N and r_regions, determine n_regions,
    #  a list of the natural number of regions in each collar and the polar caps.
    #  This list is as close as possible to r_regions, using rounding.
    #  The number of elements is n_collars+2.
    #  n_regions[1] is 1.
    #  n_regions[n_collars+2] is 1.
    #  The sum of n_regions is N.
    #
    discrepancy = 0.0
    for zone_n = 1 : (ncollars + 2)
        n_regions[zone_n] = round(Int,r_regions[zone_n] + discrepancy)
        discrepancy = discrepancy + r_regions[zone_n] - n_regions[zone_n]
    end
end

function cap_colats!(N::Int, n_collars::Int, c_polar::Float64, n_regions::Array{Int,1}, c_caps::Array{Float64,1})
    # CAP_COLATS Colatitudes of spherical caps enclosing cumulative sum of
    # regions
    #
    #  Given dim, N, c_polar and n_regions, determine c_caps,
    #  an increasing list of colatitudes of spherical caps which enclose the same
    #  area
    #  as that given by the cumulative sum of regions.
    #  The number of elements is n_collars+2.
    #  c_caps[1] is c_polar.
    #  c_caps[n_collars+1] is Pi-c_polar.
    #  c_caps[n_collars+2] is Pi.
    #
    #  c_caps = cap_colats(dim,N,c_polar,n_regions)

    c_caps[1]                = c_polar
    ideal_region_area = area_of_ideal_region(N)
    subtotal_n_regions   = 1
    for  collar_n = 1 : n_collars
        subtotal_n_regions   = subtotal_n_regions + n_regions[1 + collar_n]
        c_caps[collar_n + 1] = sradius_of_cap(subtotal_n_regions * ideal_region_area)
    end
    c_caps[n_collars + 2] = pi
end


function eq_caps(N::Int)
    #
    # eq_regions uses the zonal equal area sphere partitioning algorithm to
    # partition
    # the surface of a sphere into N regions of equal area and small diameter.
    #

    if N == 1 
        #
        # We have only one region, which must be the whole sphere.
        #
        n_regions=zeros(Int,1)
        s_cap=zeros(Float64,1)
        n_regions[1] = 1
        s_cap[1]     = pi
        # int n_regions_ns=1
    else
        #
        # Given N, determine c_polar
        # the colatitude of the North polar spherical cap.
        #
        c_polar = polar_colat(N)

        #
        # Given N, determine the ideal angle for spherical collars.
        # Based on N, this ideal angle, and c_polar,
        # determine n_collars, the number of collars between the polar caps.
        #
        n_collars = num_collars(N, c_polar, ideal_collar_angle(N))

        # int n_regions_ns=n_collars+2
        #
        # Given N, c_polar and n_collars, determine r_regions,
        # a list of the ideal real number of regions in each collar,
        # plus the polar caps.
        # The number of elements is n_collars+2.
        # r_regions[0] is 1.
        # r_regions[2+n_collars-1] is 1.
        # The sum of r_regions is N.
        r_regions=zeros(Float64,n_collars + 2)
        ideal_region_list!(N, c_polar, n_collars, r_regions)
        #
        # Given N and r_regions, determine n_regions, a list of the natural number
        # of regions in each collar and the polar caps.
        # This list is as close as possible to r_regions.
        # The number of elements is n_collars+2.
        # n_regions[0] is 1.
        # n_regions[2+n_collars-1] is 1.
        # The sum of n_regions is N.
        #
        n_regions=zeros(Int,n_collars + 2)
        round_to_naturals!(N, n_collars, r_regions, n_regions)
        #
        # Given dim, N, c_polar and n_regions, determine s_cap,
        # an increasing list of colatitudes of spherical caps which enclose the
        # same area
        # as that given by the cumulative sum of regions.
        # The number of elements is n_collars+2.
        # s_cap[0] is c_polar.
        # s_cap[n_collars]   is Pi-c_polar.
        # s_cap[n_collars+1] is Pi.
        #
        s_cap = zeros(Float64,n_collars + 2)
        cap_colats!(N, n_collars, c_polar, n_regions, s_cap)
    end
    return (n_regions,s_cap)
    # int n_regions_ew=maxval(n_regions(:))
end


function eq_regions!(N::Int, xmin::Array{Float64,1}, xmax::Array{Float64,1}, ymin::Array{Float64,1}, ymax::Array{Float64,1})
    # EQ_REGIONS Recursive zonal equal area (EQ) partition of sphere
    #
    # Syntax
    #  [regions,dim_1_rot] = eq_regions(dim,N,options)
    #
    # Description
    #  REGIONS = EQ_REGIONS(dim,N) uses the recursive zonal equal area sphere
    #  partitioning algorithm to partition S^dim (the unit sphere in dim+1
    #  dimensional space) into N regions of equal area and small diameter.
    #
    #  The arguments dim and N must be positive integers.
    #
    #  The result REGIONS is a (dim by 2 by N) array, representing the regions
    #  of S^dim. Each element represents a pair of vertex points in spherical
    #  polar
    #  coordinates.
    #
    #  Each region is defined as a product of intervals in spherical polar
    #  coordinates. The pair of vertex points regions(:,1,n) and regions(:,2,n)
    #  give
    #  the lower and upper limits of each interval.
    #
    #  REGIONS = EQ_REGIONS(dim,N,'offset','extra') uses experimental extra
    #  offsets for S^2 and S^3 to try to minimize energy. If dim > 3, extra
    #  offsets
    #  are not used.
    #
    #  REGIONS = EQ_REGIONS(dim,N,extra_offset) uses experimental extra offsets
    #  if extra_offset is true or non-zero.
    #
    #  [REGIONS,DIM_1_ROT] = EQ_REGIONS(dim,N) also returns DIM_1_ROT, a cell
    #  array containing N rotation matrices, one per region, each of size dim by
    #  dim.
    #  These describe the R^dim rotation needed to place the region in its final
    #  position.
    #
    #  [REGIONS,DIM_1_ROT] = EQ_REGIONS(dim,N,'offset','extra') partitions S^dim
    #  into N regions, using extra offsets, and also returning DIM_1_ROT, as
    #  above.
    #

    if N == 1
        #
        # We have only one region, which must be the whole sphere.
        #
        xmin[1] = 0.
        ymin[1] = -0.5 * pi
        xmax[1] = 2. * pi
        ymax[1] = 0.5 * pi
        return
    end
    #
    # Start the partition of the sphere into N regions by partitioning
    # to caps defined in the current dimension.
    #
    (n_regions,s_cap) = eq_caps(N)
    #
    # s_cap is an increasing list of colatitudes of the caps.
    #
    #
    # We have a number of zones: two polar caps and a number of collars.
    # n_regions is the list of the number of regions in each zone.
    #
    n_collars = size(n_regions,1) - 2
    #
    # Start with the top cap.
    #
    xmin[1] = 0.
    ymin[1] = 0.5 * pi - s_cap[1]
    xmax[1] = 2. * pi
    ymax[1] = 0.5 * pi

    region_n = 2
    for collar_n = 1 : n_collars
        for region_ew = 1 : n_regions[collar_n + 1]
            xmin[region_n] = 2. * pi / (n_regions[collar_n + 1]) * (region_ew - 1.0)
            ymin[region_n] = 0.5 * pi - s_cap[collar_n + 1]
            xmax[region_n] = 2. * pi / (n_regions[collar_n + 1]) * region_ew 
            ymax[region_n] = 0.5 * pi - s_cap[collar_n]
            region_n += 1
        end
    end
    #
    # End with the bottom cap.
    #
    xmin[N] = 0.
    ymin[N] = -0.5 * pi
    xmax[N] = 2. * pi
    ymax[N] = 0.5 * pi - s_cap[end - 1]
end


function compare_NS_WE(node1, node2) 
  if node1[2] > node2[2]
    return true
 end
 if node1[2] == node2[2]
   return node1[1] < node2[1]
 end
 return false
end

function compare_WE_NS(node1, node2) 
  if node1[1] < node2[1] 
    return true
  end
  if node1[1] == node2[1] 
    return node1[2] > node2[2]
  end
  return false
end

