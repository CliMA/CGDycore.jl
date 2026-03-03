function OrderConditionsRosenbrockW(RosW,Order)

# (O1) Sum bj = 1

# (O2) Sum bjajk = 1
# (O2) Sum bjdjk = 0 

# (O3) Sum bjajkajl = 1/3
# (O3) Sum bjajkakl = 1/6
# (O3) Sum bjajkdkl = 0 
# (O3) Sum bjdjkakl = 0 
# (O3) Sum bjdjkdkl = 0 

  nStage = RosW.nStage
  alpha = RosW.alpha
  gamma = RosW.gamma
  b = RosW.b
  gammaD = RosW.gammaD

  if Order == 1
    NumOrder =1
  elseif Order == 2
    NumOrder = 3
  elseif Order == 3
    NumOrder = 8
  end  
  O = zeros(NumOrder)

# O1
      
  O[1] = -1
  for i = 1 : nStage
    O[1] += b[i]
  end
  if Order == 1
    return O
  end

# O2

  O[2] = -1.0/2.0
  for i = 1 : nStage
    t = 0.0
    for j = 1 : i-1
      t += alpha[i,j]
    end
    O[2] += b[i]*t
  end

  O[3] = 0.0
  for i = 1 : nStage
    t = 0.0
    for j = 1 : i
      t += gamma[i,j]
    end
    O[3] += b[i]*t
  end
  if Order == 2
    return O
  end

# O3
# (O3) Sum bjajkajl = 1/3
  O[4] = -1.0/3.0
  for i = 1 : nStage
    t = 0
    for j = 1 : nStage
      t = t + alpha[i,j]
    end
    O[4] += b[i] * t^2
  end
# (O3) Sum bjajkakl = 1/6
  O[5] = -1.0/6.0
  for i = 1 : nStage
    t1 = 0
    for j = 1 : nStage
      t2 = 0
      for k = 1 : nStage
        t2 += alpha[j,k]
      end
      t1 += alpha[i,j] * t2
    end
    O[5] += b[i] * t1
  end

  O[6] = 0
  for i = 1 : nStage
    t1 = 0
    for j = 1 : nStage
      t2 = 0
      for k = 1 : nStage
        t2 += gamma[j,k]
      end
      t1 += alpha[i,j] * t2
    end
    O[6] +=  b[i] * t1
  end
  O[7] = 0
  for i = 1 : nStage
    t1 = 0
    for j = 1 : nStage
      t2 = 0
      for k = 1 : nStage
        t2 += alpha[j,k]
      end
      t1 += gamma[i,j] * t2
    end
    O[7] += b[i] * t1
  end
  O[8] = 0
  for i = 1 : nStage
    t1 = 0
    for j = 1 : nStage
      t2 = 0
      for k = 1 : nStage
        t2 += gamma[j,k]
      end
      t1 += gamma[i,j] * t2
    end
    O[8] += b[i] * t1
  end
  if Order == 3
    return O
  end
end




