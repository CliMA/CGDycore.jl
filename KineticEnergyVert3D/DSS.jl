function DSSC!(c,JRho)
  Nx=size(c,1)
  Ny=size(c,2)
  Nz=size(c,3)
  OPx=size(c,4)
  OPy=size(c,5)
  OPz=size(c,6)
  # Face x
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        if ix > 1
          @views cL = c[ix-1,iy,iz,OPx,:,:]
          @views JRL = JRho[ix-1,iy,iz,OPx,:,:]
        else
          @views cL = c[Nx,iy,iz,OPx,:,:]
          @views JRL = JRho[Nx,iy,iz,OPx,:,:]
        end
        @views cR = c[ix,iy,iz,1,:,:]
        @views JRR = JRho[ix,iy,iz,1,:,:]
        @inbounds for i = 1 : OPy
          @inbounds for j = 1 : OPz
            cM = cL[i,j] + cR[i,j]
            JRM =   JRL[i,j] + JRR[i,j]
            cL[i,j] = cM
            cR[i,j] = cM
            JRL[i,j] = JRM
            JRR[i,j] = JRM
          end
        end
      end
    end
  end
# y
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        if iy > 1
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz
              cM = c[ix,iy-1,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              JRM = JRho[ix,iy-1,iz,i,OPy,j] + JRho[ix,iy,iz,i,1,j]
              c[ix,iy-1,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              JRho[ix,iy-1,iz,i,OPy,j] = JRM
              JRho[ix,iy,iz,i,1,j] = JRM
            end
          end
        else
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz
              cM = c[ix,Ny,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              JRM = JRho[ix,Ny,iz,i,OPy,j] + JRho[ix,iy,iz,i,1,j]
              c[ix,Ny,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              JRho[ix,Ny,iz,i,OPy,j] = JRM
              JRho[ix,iy,iz,i,1,j] = JRM
            end
          end
        end
      end
    end
  end

  @. c /= JRho
end

function AverageC!(c)
  Nx=size(c,1)
  Ny=size(c,2)
  Nz=size(c,3)
  OPx=size(c,4)
  OPy=size(c,5)
  OPz=size(c,6)
  Mass = similar(c)
  @. Mass = 1.0
  # Face x
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        if ix > 1
          @views cL = c[ix-1,iy,iz,OPx,:,:]
          @views MassL = Mass[ix-1,iy,iz,OPx,:,:]
        else
          @views cL = c[Nx,iy,iz,OPx,:,:]
          @views MassL = Mass[Nx,iy,iz,OPx,:,:]
        end
        @views cR = c[ix,iy,iz,1,:,:]
        @views MassR = Mass[ix,iy,iz,1,:,:]
        @inbounds for i = 1 : OPy
          @inbounds for j = 1 : OPz
            cM = cL[i,j] + cR[i,j]
            MassM =   MassL[i,j] + MassR[i,j]
            cL[i,j] = cM
            cR[i,j] = cM
            MassL[i,j] = MassM
            MassR[i,j] = MassM
          end
        end
      end
    end
  end
# y
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        if iy > 1
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz
              cM = c[ix,iy-1,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              MassM = Mass[ix,iy-1,iz,i,OPy,j] + Mass[ix,iy,iz,i,1,j]
              c[ix,iy-1,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              Mass[ix,iy-1,iz,i,OPy,j] = MassM
              Mass[ix,iy,iz,i,1,j] = MassM
            end
          end
        else
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz
              cM = c[ix,Ny,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              MassM = Mass[ix,Ny,iz,i,OPy,j] + Mass[ix,iy,iz,i,1,j]
              c[ix,Ny,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              Mass[ix,Ny,iz,i,OPy,j] = MassM
              Mass[ix,iy,iz,i,1,j] = MassM
            end
          end
        end
      end
    end
  end

  @. c /= Mass
end

function DSSF!(c,JRho)
  Nx=size(c,1)
  Ny=size(c,2)
  Nz=size(c,3)
  OPx=size(c,4)
  OPy=size(c,5)
  OPz=size(c,6)
  # Face x
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx  
        if ix > 1
          @views cL = c[ix-1,iy,iz,OPx,:,:]  
          @views JRL = JRho[ix-1,iy,iz,OPx,:,:]  
        else
          @views cL = c[Nx,iy,iz,OPx,:,:]  
          @views JRL = JRho[Nx,iy,iz,OPx,:,:]  
        end  
        @views cR = c[ix,iy,iz,1,:,:]  
        @views JRR = JRho[ix,iy,iz,1,:,:]  
        @inbounds for i = 1 : OPy  
          @inbounds for j = 1 : OPz   
            cM = cL[i,j] + cR[i,j]
            JRM =   JRL[i,j] + JRR[i,j] 
            cL[i,j] = cM
            cR[i,j] = cM
            JRL[i,j] = JRM
            JRR[i,j] = JRM
          end
        end
      end
    end
  end
  # y
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        if iy > 1
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz  
              cM = c[ix,iy-1,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              JRM = JRho[ix,iy-1,iz,i,OPy,j] + JRho[ix,iy,iz,i,1,j]
              c[ix,iy-1,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              JRho[ix,iy-1,iz,i,OPy,j] = JRM
              JRho[ix,iy,iz,i,1,j] = JRM
            end
          end  
        else
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz  
              cM = c[ix,Ny,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              JRM = JRho[ix,Ny,iz,i,OPy,j] + JRho[ix,iy,iz,i,1,j]
              c[ix,Ny,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              JRho[ix,Ny,iz,i,OPy,j] = JRM
              JRho[ix,iy,iz,i,1,j] = JRM
            end
          end  
        end
      end
    end
  end
  # z
  @inbounds for iz = 2 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        @views cL = c[ix,iy,iz-1,:,:,OPz]
        @views JRL = JRho[ix,iy,iz-1,:,:,OPz]
        @views cR = c[ix,iy,iz,:,:,1]
        @views JRR = JRho[ix,iy,iz,:,:,1]
        @inbounds for i = 1 : OPx
          @inbounds for j = 1 : OPy
            cM = cL[i,j] + cR[i,j]
            JRM =   JRL[i,j] + JRR[i,j]
            cL[i,j] = cM
            cR[i,j] = cM
            JRL[i,j] = JRM
            JRR[i,j] = JRM
          end
        end
      end
    end
  end
  @. c /= JRho

end  

function AverageF!(c)
  Nx=size(c,1)
  Ny=size(c,2)
  Nz=size(c,3)
  OPx=size(c,4)
  OPy=size(c,5)
  OPz=size(c,6)
  Mass = similar(c)
  @. Mass = 1.0
  # Face x
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx  
        if ix > 1
          @views cL = c[ix-1,iy,iz,OPx,:,:]  
          @views MassL = Mass[ix-1,iy,iz,OPx,:,:]  
        else
          @views cL = c[Nx,iy,iz,OPx,:,:]  
          @views MassL = Mass[Nx,iy,iz,OPx,:,:]  
        end  
        @views cR = c[ix,iy,iz,1,:,:]  
        @views MassR = Mass[ix,iy,iz,1,:,:]  
        @inbounds for i = 1 : OPy  
          @inbounds for j = 1 : OPz   
            cM = cL[i,j] + cR[i,j]
            MassM =   MassL[i,j] + MassR[i,j] 
            cL[i,j] = cM
            cR[i,j] = cM
            MassL[i,j] = MassM
            MassR[i,j] = MassM
          end
        end
      end
    end
  end
  # y
  @inbounds for iz = 1 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        if iy > 1
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz  
              cM = c[ix,iy-1,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              MassM = Mass[ix,iy-1,iz,i,OPy,j] + Mass[ix,iy,iz,i,1,j]
              c[ix,iy-1,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              Mass[ix,iy-1,iz,i,OPy,j] = MassM
              Mass[ix,iy,iz,i,1,j] = MassM
            end
          end  
        else
          @inbounds for i = 1 : OPx
            @inbounds for j = 1 : OPz  
              cM = c[ix,Ny,iz,i,OPy,j] + c[ix,iy,iz,i,1,j]
              MassM = Mass[ix,Ny,iz,i,OPy,j] + Mass[ix,iy,iz,i,1,j]
              c[ix,Ny,iz,i,OPy,j] = cM
              c[ix,iy,iz,i,1,j] = cM
              Mass[ix,Ny,iz,i,OPy,j] = MassM
              Mass[ix,iy,iz,i,1,j] = MassM
            end
          end  
        end
      end
    end
  end
  # z
  @inbounds for iz = 2 : Nz
    @inbounds for iy = 1 : Ny
      @inbounds for ix = 1 : Nx
        @views cL = c[ix,iy,iz-1,:,:,OPz]
        @views MassL = Mass[ix,iy,iz-1,:,:,OPz]
        @views cR = c[ix,iy,iz,:,:,1]
        @views MassR = Mass[ix,iy,iz,:,:,1]
        @inbounds for i = 1 : OPx
          @inbounds for j = 1 : OPy
            cM = cL[i,j] + cR[i,j]
            MassM =   MassL[i,j] + MassR[i,j]
            cL[i,j] = cM
            cR[i,j] = cM
            MassL[i,j] = MassM
            MassR[i,j] = MassM
          end
        end
      end
    end
  end
  @. c /= Mass

end  

function AverageFXY!(c,J)
  Nx=size(c,1)
  Ny=size(c,2)
  OPx=size(c,3)
  OPy=size(c,4)
  # Face x
  @inbounds for iy = 1 : Ny
    @inbounds for ix = 1 : Nx  
      if ix > 1
        @views cL = c[ix-1,iy,OPx,:]  
        @views JL = J[ix-1,iy,OPx,:]  
      else
        @views cL = c[Nx,iy,OPx,:]  
        @views JL = J[Nx,iy,OPx,:]  
      end  
      @views cR = c[ix,iy,1,:]  
      @views JR = J[ix,iy,1,:]  
      @inbounds for i = 1 : OPy  
        cM = cL[i] + cR[i]
        JM = JL[i] + JR[i]
        cL[i] = cM
        cR[i] = cM
        JL[i] = JM
        JR[i] = JM
      end
    end
  end
  # y
  @inbounds for iy = 1 : Ny
    @inbounds for ix = 1 : Nx
      if iy > 1
        @inbounds for i = 1 : OPx  
          cM = c[ix,iy-1,i,OPy] + c[ix,iy,i,1]
          JM = J[ix,iy-1,i,OPy] + J[ix,iy,i,1]
          c[ix,iy-1,i,OPy] = cM
          c[ix,iy,i,1] = cM
          J[ix,iy-1,i,OPy] = JM
          J[ix,iy,i,1] = JM
        end  
      else
        @inbounds for i = 1 : OPx  
          cM = c[ix,Ny,i,OPy] + c[ix,iy,i,1]
          JM = J[ix,Ny,i,OPy] + J[ix,iy,i,1]
          c[ix,Ny,i,OPy] = cM
          c[ix,iy,i,1] = cM
          J[ix,Ny,i,OPy] = JM
          J[ix,iy,i,1] = JM
        end  
      end
    end
  end
  @. c /= J
end  


