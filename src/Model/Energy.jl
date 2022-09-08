function Energy!(E,RhoTh,Rho,Tr,V,CG,Global)
  (; Rd,
     Cvd,
     Cpd,
     Rv,
     Cvv,
     Cpv,
     Cpl,
     p0,
     Grav,
     kappa) = Global.Phys


  @views RhoCG = Global.Cache.RhoCG[:,:,:]
  @views ThCG = Global.Cache.ThCG[:,:,:]
  @views ECG = Global.Cache.ThCG[:,:,:]
  v1CG = Global.Cache.v1CG
  v2CG = Global.Cache.v2CG
  wCG = Global.Cache.wCG
  wCCG = Global.Cache.wCCG
  KE = Global.Cache.KE
  JC = Global.Metric.JC
  Equation = Global.Model.Equation
  OP = CG.OrdPoly + 1
  nz = Global.Grid.nz
  NF = Global.Grid.NumFaces
  NumTr = Global.Model.NumTr
  zP = Global.Metric.zP
  XP = zeros(3)
  @show zP[:,1]
  if Global.Model.Equation == "Compressible"
    if Global.Model.Thermo == "TotalEnergy"  
      @inbounds for iF = 1 : NF
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
            @inbounds for iz=1:nz
              RhoCG[iP,jP,iz] = Rho[iz,ind]
              v1CG[iP,jP,iz] = V[iz,ind,1]
              v2CG[iP,jP,iz] = V[iz,ind,2]
              wCG[iP,jP,iz+1] = V[iz,ind,3]
              ThCG[iP,jP,iz] = RhoTh[iz,ind]
              @inbounds for iT = 1:NumTr
                TrCG[iP,jP,iz,iT] = Tr[iz,ind,iT]
              end
            end
          end
        end
#       Diagnostic values
#       Boundary value for vertical velocity and cell center
        BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
#       Kinetic energy
        @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
#       Pressure
        @. ECG = p0 * (Rd * ThCG / p0)^(1.0 / (1.0 - kappa))
        @. ECG = ECG / (Rd *  RhoCG)
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
            @inbounds for iz=1:nz
              E[iz,ind] += RhoCG[iP,jP,iz] * (Cvd * ECG[iP,jP,iz] + KE[iP,jP,iz] + Grav*zP[iz,ind]) * JC[iP,jP,iz,iF] / CG.M[iz,ind] 
            end
          end
        end
      end
    elseif Global.Model.Thermo == "InternalEnergy"  
      @inbounds for iF = 1 : NF
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
            @inbounds for iz=1:nz
              RhoCG[iP,jP,iz] = Rho[iz,ind]
              ThCG[iP,jP,iz] = RhoTh[iz,ind]
              @inbounds for iT = 1:NumTr
                TrCG[iP,jP,iz,iT] = Tr[iz,ind,iT]
              end
            end
          end
        end
#       Diagnostic values
#       Pressure
        @inbounds for jP=1:OP
          @inbounds for iP=1:OP
            ind = CG.Glob[iP,jP,iF]
            @inbounds for iz=1:nz
              EE = p0 * (Rd * RhoTh[iz,ind] / p0)^(1.0 / (1.0 - kappa))
              E[iz,ind] += Cvd / Rd * EE  * JC[iP,jP,iz,iF] / CG.M[iz,ind] 
            end
          end
        end
      end
    end
    ExchangeData!(E,Global.Exchange)
  end  
end  

      
