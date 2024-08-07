function QuadTriangle(Order)
  if Order == 2
    w = zeros(3)
    ksi = zeros(3,2)
    w .= 0.333333333333333
    ksi[1,1] = 0.166666666666667
    ksi[1,2] = 0.666666666666667
    ksi[2,1] = 0.166666666666667
    ksi[2,2] = 0.166666666666667
    ksi[3,1] = 0.666666666666667
    ksi[3,2] = 0.166666666666667
    @. w *= 2
    @. ksi = 2 * ksi - 1 
  elseif Order == 4
    w = zeros(6)
    ksi = zeros(6,2)
    w[1] = 0.223381589678011
    w[2] = 0.223381589678011
    w[3] = 0.223381589678011
    w[4] = 0.109951743655322
    w[5] = 0.109951743655322
    w[6] = 0.109951743655322

    ksi[1,1] = 0.445948490915965
    ksi[2,1] = 0.445948490915965
    ksi[3,1] = 0.108103018168070
    ksi[4,1] = 0.091576213509771
    ksi[5,1] = 0.091576213509771
    ksi[6,1] = 0.816847572980459

    ksi[1,2] = .108103018168070
    ksi[2,2] = 0.445948490915965
    ksi[3,2] = 0.445948490915965
    ksi[4,2] = 0.816847572980459
    ksi[5,2] = 0.091576213509771
    ksi[6,2] = 0.091576213509771
    @. w *= 2
    @. ksi = 2 * ksi - 1 
  end
  return w, ksi
end  
