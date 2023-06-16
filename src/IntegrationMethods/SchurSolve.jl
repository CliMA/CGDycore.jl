function triSolve!(x,tri,b)
# Thomas algorithm, b and x used as intermediate storage
# b is destroyed after the elemination
  n = size(tri,2)
  x[1] = tri[2,1]
  @inbounds for i=1:n-1
    @fastmath x[i+1] = tri[2,i+1]-tri[3,i]/x[i]*tri[1,i+1]  
    @fastmath b[i+1] = b[i+1]-tri[3,i]/x[i]*b[i]  
  end
  @fastmath x[n] = b[n]/x[n]
  @inbounds for i=n-1:-1:1
    @fastmath x[i]=(b[i] - tri[1,i+1]*x[i+1])/x[i]  
  end
end

function mulUL!(tri,biU,biL)
  n = size(biL,2)
  @fastmath tri[2,1] = tri[2,1] - biU[2,1]*biL[1,1] - biU[1,2]*biL[2,1]
  @fastmath tri[3,1] = tri[3,1] - biU[2,2]*biL[2,1]
  @inbounds for i=2:n-1
    @fastmath tri[1,i] = tri[1,i] - biU[1,i]*biL[1,i]
    @fastmath tri[2,i] = tri[2,i] - biU[2,i]*biL[1,i] - biU[1,i+1]*biL[2,i]
    @fastmath tri[3,i] = tri[3,i] - biU[2,i+1]*biL[2,i]
  end
  i = n
  @fastmath tri[1,i] = tri[1,i] - biU[1,i]*biL[1,i]
  @fastmath tri[2,i] = tri[2,i] - biU[2,i]*biL[1,i]
end

function mulbiUv!(u,biU,v)
  n = size(biU,2)
  @inbounds for i=1:n-1
    @fastmath u[i] = u[i] + biU[2,i]*v[i] + biU[1,i+1]*v[i+1]
  end
  @fastmath u[n] = u[n] + biU[2,n]*v[n]
end

function mulbiLv!(u,biL,v)
  n = size(biL,2)
  @fastmath u[1] = u[1] + biL[1,1]*v[1]
  @inbounds for i=2:n
    @fastmath u[i] = u[i] + biL[2,i-1]*v[i-1] + biL[1,i]*v[i]
  end
end

function SchurSolve!(k,v,J,fac,Global)
#   sw=(spdiags(repmat(invfac2,n,1),0,n,n)-invfac*JWW-JWRho*JRhoW-JWTh*JThW)\
#     (invfac*rw+JWRho*rRho+JWTh*rTh)
  n1 = size(v,1)
  n2 = size(v,2)
  n = n1 * n2
  tri = J.tri
  JWRho = J.JWRho
  JRhoW = J.JRhoW
  JWTh = J.JWTh
  JThW=J.JThW
  JTrW=J.JTrW
  JWW=J.JWW
  JDiff = J.JDiff
  JAdvC = J.JAdvC
  JAdvF = J.JAdvF
  @views CdTh = Global.Cache.Aux2DG[:,:,1]
  @views CdTr = Global.Cache.Aux2DG[:,:,2:end]
  NumV = Global.Model.NumV
  if size(k,3) > Global.Model.NumV
    NumTr = Global.Model.NumTr
  else
    NumTr = 0
  end  

  invfac=1/fac
  invfac2=invfac/fac

  if Global.Model.VerticalDiffusion
    @inbounds for in2=1:n2
      if J.CompTri
        @views @. JAdvC[2,:,in2] = invfac  + JAdvC[2,:,in2]  
        @views @. JAdvF[2,:,in2] = invfac  + JAdvF[2,:,in2]  
        @views @. JDiff[:,:,in2] += JAdvC[:,:,in2]  
      end
      JDiff[2,1,in2] += CdTh[1,in2]  
      @views rTh = v[:,in2,5]
      @views sTh = k[:,in2,5]
      @views triSolve!(sTh,JDiff[:,:,in2],rTh)
      @. rTh = invfac * sTh
      JDiff[2,1,in2] -= CdTh[1,in2] 
      @views rTh = v[:,in2,2]
      @views sTh = k[:,in2,2]
      @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
      @. rTh = invfac * sTh
      @views rTh=v[:,in2,3]
      @views sTh=k[:,in2,3]
      @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
      @. rTh = invfac * sTh
      @views rTh=v[2:n1,in2,4]
      @views sTh=k[2:n1,in2,4]
      @views triSolve!(sTh,JAdvF[:,:,in2],rTh)
      @. rTh = invfac * sTh
    end
    if Global.Model.Equation == "CompressibleMoist"
      @inbounds for in2=1:n2
        JDiff[2,1,in2] += CdTr[1,in2,1] 
        @views rTh=v[:,in2,NumV+1]
        @views sTh=k[:,in2,NumV+1]
        @views triSolve!(sTh,JDiff[:,:,in2],rTh)
        @. rTh = invfac * sTh
        JDiff[2,1,in2] -= CdTr[1,in2,1] 
      end  
    end    
  end    

  if Global.Model.Equation == "Compressible"
    @inbounds for in2=1:n2
      @views rRho=v[:,in2,1]
      @views rTh=v[:,in2,5]
      @views rw=v[:,in2,4]
      @views sw=k[:,in2,4]
      if Global.Model.Damping
        if J.CompTri
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2  - invfac * JWW[1,:,in2]
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWTh[:,:,in2],JThW[:,:,in2])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      else
        if J.CompTri  
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2   
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWTh[:,:,in2],JThW[:,:,in2])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      end
      @views mulbiLv!(rRho,JRhoW[:,:,in2],sw)
      @views mulbiLv!(rTh,JThW[:,:,in2],sw)

      @views @. k[:,in2,1] = fac * rRho
      @views @. k[:,in2,2:3] = fac * v[:,in2,2:3]
      @views @. k[:,in2,5] = fac * rTh
      @inbounds for iT = 1 : NumTr
        @views mulbiLv!(v[:,in2,5+iT],JTrW[:,:,in2,iT],sw)  
        @views @. k[:,in2,5+iT] = fac * v[:,in2,5+iT]
      end  
      if Global.Model.Damping
        @views @. sw = sw / (1.0 - invfac * JWW[1,:,in2])   
      end  
    end
  elseif Global.Model.Equation == "CompressibleMoist"
    RhoVPos = Global.Model.RhoVPos
    JWRhoV=J.JWRhoV
    @inbounds for in2=1:n2
      @views rRho=v[:,in2,1]
      @views rTh=v[:,in2,5]
      @views rRhoV=v[:,in2,NumV + RhoVPos]
      @views rw=v[:,in2,4]
      @views sw=k[:,in2,4]
      invfac=1/fac
      invfac2=invfac/fac
      if Global.Model.Damping
        if J.CompTri
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2  - invfac * JWW[1,:,in2]
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWTh[:,:,in2],JThW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoV[:,:,in2],JTrW[:,:,in2,RhoVPos])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWRhoV[:,:,in2],rRhoV)
        @views mulbiUv!(rw,JWTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      else
        if J.CompTri  
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2   
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWTh[:,:,in2],JThW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoV[:,:,in2],JTrW[:,:,in2,RhoVPos])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWRhoV[:,:,in2],rRhoV)
        @views mulbiUv!(rw,JWTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      end
      @views mulbiLv!(rRho,JRhoW[:,:,in2],sw)
      @views mulbiLv!(rTh,JThW[:,:,in2],sw)

      @views @. k[:,in2,1] = fac * rRho
      @views @. k[:,in2,2:3] = fac * v[:,in2,2:3]
      @views @. k[:,in2,5] = fac * rTh
      @inbounds for iT = 1 : NumTr
        @views mulbiLv!(v[:,in2,5+iT],JTrW[:,:,in2,iT],sw)  
        @views @. k[:,in2,5+iT] = fac * v[:,in2,5+iT]
      end  
      if Global.Model.Damping
        @views @. sw = sw / (1.0 - invfac * JWW[1,:,in2])   
      end  
    end
  end
  J.CompTri=false
end


