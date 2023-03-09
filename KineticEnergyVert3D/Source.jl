function Damping!(FwF,wF,X,Fe)
  H = 15600.0
  StrideDamp = 10000.0
  Relax = 0.0 # 0.1
  Nz = size(FwF,1)
  OrdPolyX = Fe.OrdPolyX
  OrdPolyY = Fe.OrdPolyY
  OrdPolyZ = Fe.OrdPolyZ
  @inbounds for iz = 1 : Nz
    @inbounds for i = 1 : OrdPolyX + 1
      @inbounds for j = 1 : OrdPolyY + 1
        @inbounds for k = 1 : OrdPolyZ + 1
          zLoc = X[iz,i,j,k,3]
          if zLoc >= H - StrideDamp
            Damp = Relax * sin(0.5*pi*(1.0 - (H - zLoc)/StrideDamp))^2;
            FwF[iz,i,j,k] -= Damp * wF[iz,i,j,k]
          end
        end
      end
    end
  end
end  

function Buoyancy!(FwF,RhoF,J,Phys)

  @. FwF -= Phys.Grav * RhoF * J

end
