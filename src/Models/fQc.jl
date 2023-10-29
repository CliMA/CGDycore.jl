function fQc(x,time,Global,Param,Profile)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.Problem)
  if str == "bryanfritschcart"
    z = x[3]
    @views zP = Profile[:,1]
    iz = 1000
    for i = 2:size(zP,1)
      if z <= zP[i]
        iz = i - 1 
        break
      end  
    end 
    z_l = zP[iz]
    Rho_l = Profile[iz,2]
    Theta_l = Profile[iz,3]
    RhoV_l = Profile[iz,4]
    RhoC_l = Profile[iz,5]
    z_r = zP[iz+1]
    Rho_r = Profile[iz+1,2]
    Theta_r = Profile[iz+1,3]
    RhoV_r = Profile[iz+1,4]
    RhoC_r = Profile[iz+1,5]
    Rho = (Rho_r * (z - z_l) + Rho_l * (z_r - z)) / (z_r - z_l)
    Theta = (Theta_r * (z - z_l) + Theta_l * (z_r - z)) / (z_r - z_l)
    RhoV = (RhoV_r * (z - z_l) + RhoV_l * (z_r - z)) / (z_r - z_l)
    RhoC = (RhoC_r * (z - z_l) + RhoC_l * (z_r - z)) / (z_r - z_l)

    Rho, Theta, qv, qc = PerturbMoistProfile(x, Rho, Rho*Theta, RhoV, RhoC, Phys, Param)
  else  
    qc = 0  
  end
  return qc
end    
