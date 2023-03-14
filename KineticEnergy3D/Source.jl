function Damping!(FwF,wF,X,Fe)
  H = 15600.0
  StrideDamp = 10000.0
  Relax = 0.1
  Nz = size(FwF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  @inbounds for iz = 1 : Nz
    @inbounds for i = 1 : OrdPolyX + 1
      @inbounds for j = 1 : OrdPolyY + 1
        if iz == 1
          zLoc = X[iz,i,j,1,3]
        else
          zLoc = X[iz-1,i,j,2,3]  
        end  
        if zLoc >= H - StrideDamp
          Damp = Relax * sin(0.5*pi*(1.0 - (H - zLoc)/StrideDamp))^2;
          FwF[iz,i,j] -= Damp * wF[iz,i,j]
        end
      end
    end
  end
end  

function Buoyancy!(FwF,RhoC,J,Phys)

  @views @. FwF[2:end-1,:,:] -= Phys.Grav * 
    (RhoC[1:end-1,:,:] * J[1:end-1,:,:,2] +
    RhoC[2:end,:,:] * J[2:end,:,:,1])

end
