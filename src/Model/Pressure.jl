function Pressure(U,KE,zP,Global)
  Rd = Global.Phys.Rd
  Cvd = Global.Phys.Cvd
  Grav = Global.Phys.Grav
  if Global.Model.Equation == "Compressible"
    p=(Rd/Cvd)*(U[5]-U[1]*(KE+Grav*zP))
  elseif Global.Model.Equation == "CompressibleMoist"
  end
end 

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
fast_pow(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))

function Pressure!(p::AbstractArray{Float64,3},RhoTh::AbstractArray{Float64,3},Rho::AbstractArray{Float64,3},
  Tr::AbstractArray{Float64,4},KE::AbstractArray{Float64,3},zP::AbstractArray{Float64,3},Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     Grav,
     L00,
     kappa) = Global.Phys

  
  Equation = Global.Model.Equation
  iE1 = size(p,1)
  iE2 = size(p,2)
  iE3 = size(p,3)
  if Equation == "Compressible"
     if Global.Model.Thermo == "TotalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * (RhoTh[i] - Rho[i] * (KE[i] + Grav * zP[i]))
       end  
     elseif Global.Model.Thermo == "InternalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * RhoTh[i] 
       end  
     else
       @inbounds for i in eachindex(p)  
         p[i] = p0 * fast_pow(Rd * RhoTh[i] / p0, 1.0 / (1.0 - kappa));
       end  
    end
  elseif Equation == "CompressibleMoist"
    @views TrRhoV = Tr[:,:,:,Global.Model.RhoVPos]
    @views TrRhoC = Tr[:,:,:,Global.Model.RhoCPos]
    if Global.Model.Thermo == "TotalEnergy"
    elseif Global.Model.Thermo == "InternalEnergy"
      @inbounds for i3 = 1 : iE3
        @inbounds for i2 = 1 : iE2
          @inbounds for i1 = 1 : iE1
            RhoV = TrRhoV[i1,i2,i3]
            RhoC = TrRhoC[i1,i2,i3]
            RhoD = Rho[i] - RhoV - RhoC
            p[i1,i2,i3] = (Rd * RhoD + Rv * RhoV)/(Cvd * RhoD + Cvv * RhoV + Cpl * RhoC) *
              (RhoV - L00 * RhoC)
          end    
        end    
      end    
    else
      @inbounds for i3 = 1 : iE3
        @inbounds for i2 = 1 : iE2
          @inbounds for i1 = 1 : iE1
            RhoV = TrRhoV[i1,i2,i3]
            RhoC = TrRhoC[i1,i2,i3]
            RhoD = Rho[i1,i2,i3] - RhoV - RhoC
            Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
            Rm  = Rd * RhoD + Rv * RhoV
            kappaM = Rm / Cpml
            p[i1,i2,i3] = (Rd * RhoTh[i1,i2,i3] / p0^kappaM)^(1.0 / (1.0 - kappaM))
          end  
        end  
      end  
    end  
  elseif Equation == "Shallow"
    @inbounds for i in eachindex(p)  
      p[i] = 0.5 * Grav * RhoTh[i]^2;
    end  
  end
end

function Pressure!(p::AbstractArray{Float64,2},RhoTh::AbstractArray{Float64,2},Rho::AbstractArray{Float64,2},
  Tr::AbstractArray{Float64,3},KE::AbstractArray{Float64,2},zP::AbstractArray{Float64,2},Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     Grav,
     L00,
     kappa) = Global.Phys

  
  Equation = Global.Model.Equation
  iE1 = size(p,1)
  iE2 = size(p,2)
  if Equation == "Compressible"
     if Global.Model.Thermo == "TotalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * (RhoTh[i] - Rho[i] * (KE[i] + Grav * zP[i]))
       end  
     elseif Global.Model.Thermo == "InternalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * RhoTh[i] 
       end  
     else
       @inbounds for i in eachindex(p)  
         p[i] = p0 * fast_pow(Rd * RhoTh[i] / p0, 1.0 / (1.0 - kappa));
       end  
    end
  elseif Equation == "CompressibleMoist"
    @views TrRhoV = Tr[:,:,Global.Model.RhoVPos]
    @views TrRhoC = Tr[:,:,Global.Model.RhoCPos]
    if Global.Model.Thermo == "TotalEnergy"
    elseif Global.Model.Thermo == "InternalEnergy"
      @inbounds for i2 = 1 : iE2
        @inbounds for i1 = 1 : iE1
          RhoV = TrRhoV[i1,i2]
          RhoC = TrRhoC[i1,i2]
          RhoD = Rho[i] - RhoV - RhoC
          p[i1,i2] = (Rd * RhoD + Rv * RhoV)/(Cvd * RhoD + Cvv * RhoV + Cpl * RhoC) *
            (RhoV - L00 * RhoC)
        end    
      end    
    else
      @inbounds for i2 = 1 : iE2
        @inbounds for i1 = 1 : iE1
          RhoV = TrRhoV[i1,i2]
          RhoC = TrRhoC[i1,i2]
          RhoD = Rho[i1,i2] - RhoV - RhoC
          Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
          Rm  = Rd * RhoD + Rv * RhoV
          kappaM = Rm / Cpml
          p[i1,i2] = (Rd * RhoTh[i1,i2] / p0^kappaM)^(1.0 / (1.0 - kappaM))
        end  
      end  
    end  
  elseif Equation == "Shallow"
    @inbounds for i in eachindex(p)  
      p[i] = 0.5 * Grav * RhoTh[i]^2;
    end  
  end
end

function Pressure!(p::AbstractArray{Float64,1},RhoTh::AbstractArray{Float64,1},Rho::AbstractArray{Float64,1},
  Tr::AbstractArray{Float64,2},KE::AbstractArray{Float64,1},zP::AbstractArray{Float64,1},Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     Grav,
     L00,
     kappa) = Global.Phys

  
  Equation = Global.Model.Equation
  iE1 = size(p,1)
  if Equation == "Compressible"
     if Global.Model.Thermo == "TotalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * (RhoTh[i] - Rho[i] * (KE[i] + Grav * zP[i]))
       end  
     elseif Global.Model.Thermo == "InternalEnergy"
       @inbounds for i in eachindex(p)  
         p[i] = (Rd / Cvd) * RhoTh[i] 
       end  
     else
       @inbounds for i in eachindex(p)  
         p[i] = p0 * fast_pow(Rd * RhoTh[i] / p0, 1.0 / (1.0 - kappa));
       end  
    end
  elseif Equation == "CompressibleMoist"
    @views TrRhoV = Tr[:,Global.Model.RhoVPos]
    @views TrRhoC = Tr[:,Global.Model.RhoCPos]
    if Global.Model.Thermo == "TotalEnergy"
    elseif Global.Model.Thermo == "InternalEnergy"
      @inbounds for i1 = 1 : iE1
        RhoV = TrRhoV[i1]
        RhoC = TrRhoC[i1]
        RhoD = Rho[i1] - RhoV - RhoC
        p[i1,i2] = (Rd * RhoD + Rv * RhoV)/(Cvd * RhoD + Cvv * RhoV + Cpl * RhoC) *
          (RhoV - L00 * RhoC)
      end    
    else
      @inbounds for i1 = 1 : iE1
        RhoV = TrRhoV[i1]
        RhoC = TrRhoC[i1]
        RhoD = Rho[i1] - RhoV - RhoC
        Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
        Rm  = Rd * RhoD + Rv * RhoV
        kappaM = Rm / Cpml
        if RhoTh[i1] < 0.0
          @show i1,RhoTh[i1]
          @show zP
        end  
        p[i1] = (Rd * RhoTh[i1] / p0^kappaM)^(1.0 / (1.0 - kappaM))
      end  
    end  
  elseif Equation == "Shallow"
    @inbounds for i in eachindex(p)  
      p[i] = 0.5 * Grav * RhoTh[i]^2;
    end  
  end
end

function Temperature!(T,RhoTh,Rho,Tr,Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     kappa) = Global.Phys

  Equation = Global.Model.Equation
  if Equation == "Compressible"
    if Global.Model.Thermo == "Energy"
      T=(Rd/Cvd)*(RhoTh-Rho.*(KE+Grav*repmat(Grid.zP,1,size(Rho,1))'));
    else
      @. T = p0 * (Rd * RhoTh / p0)^(1.0e0 / (1.0e0-kappa));
    end
  elseif Equation == "CompressibleMoist"
    @views @. T = PressureMoist(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
      Tr[:,Global.Model.RhoCPos],Rd,Cpd,Rv,Cpv,Cpl,p0) / (Rho * Rd + Tr[:,Global.Model.RhoVPos] * Rv)
  end
end

function PressureMoist(RhoTh,Rho,RhoV,RhoC,Rd,Cpd,Rv,Cpv,Cpl,p0)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm  = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  p = (Rd * RhoTh / p0^kappaM)^(1.0 / (1.0 - kappaM))
end


function dPresdTh!(dpdTh,RhoTh,Rho,Tr,Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     kappa) = Global.Phys

  Equation = Global.Model.Equation
  if Equation == "Compressible"
    if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
      @. dpdTh = Rd / Cvd  
    else  
      @. dpdTh=Rd*(Rd*RhoTh/p0)^(kappa/(1.0-kappa));
    end  
  elseif Equation == "CompressibleMoist"
    if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
      @. dpdTh = Rd / Cvd  
    else  
      @views @. dpdTh = dPressureMoistdTh(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
        Tr[:,Global.Model.RhoCPos],Rd,Cpd,Rv,Cpv,Cpl,p0)
    end    
  end  
end


function dPresdRhoV!(dpdRhoV,RhoTh,Rho,Tr,Pres,Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     kappa) = Global.Phys

  Equation = Global.Model.Equation
  if Equation == "Compressible"
    @. dpdRhoV = 0.0  
  elseif Equation == "CompressibleMoist"
    if Global.Model.Thermo == "TotalEnergy" || Global.Model.Thermo == "InternalEnergy"
      @. dpdRhoV = Rd / Cvd  
    else  
      @views @. dpdRhoV = dPressureMoistdRhoV(RhoTh,Rho,Tr[:,Global.Model.RhoVPos],
        Tr[:,Global.Model.RhoCPos],Pres,Rd,Cpd,Rv,Cpv,Cpl,p0)
    end    
  end  
end

function dPressureMoistdTh(RhoTh,Rho,RhoV,RhoC,Rd,Cpd,Rv,Cpv,Cpl,p0)
  RhoD = Rho - RhoV - RhoC
  Cpml = Cpd * RhoD + Cpv * RhoV + Cpl * RhoC
  Rm  = Rd * RhoD + Rv * RhoV
  kappaM = Rm / Cpml
  dpdTh=Rd*(Rd*RhoTh/p0)^(kappaM/(1-kappaM));
end

function dPressureMoistdRhoV(RhoTh,Rho,RhoV,RhoC,Pres,Rd,Cpd,Rv,Cpv,Cpl,p0)
  RhoD = Rho - RhoV - RhoC
  Rm  = Rd * RhoD + Rv * RhoV
  Cpml = RhoD * Cpd + RhoV * Cpv + RhoC * Cpl
  kappaM = Rm / Cpml
  dpdRhoV=Pres / (1.0 - kappaM)^2 * log(Rd * RhoTh / p0) *
    ((Rv - Rd) - kappaM * (Cpv - Cpd)) / Cpml 

end


