@inline function triSolve!(x,tri,b)
# Thomas algorithm, b and x used as intermediate storage
# b is destroyed after the elemination
  n = size(tri,2)
  x[1] = tri[2,1]
  for i=1:n-1
    x[i+1] = tri[2,i+1]-tri[3,i]/x[i]*tri[1,i+1]  
    b[i+1] = b[i+1]-tri[3,i]/x[i]*b[i]  
  end
  x[n] = b[n]/x[n]
  for i=n-1:-1:1
    x[i]=(b[i] - tri[1,i+1]*x[i+1])/x[i]  
  end
end

@inline function mulUL!(tri,biU,biL)
  n = size(biL,2)
  tri[2,1] = tri[2,1] - biU[2,1]*biL[1,1] - biU[1,1]*biL[2,1]
  tri[3,1] = tri[3,1] - biU[2,2]*biL[2,1]
  for i=2:n-1
    tri[1,i] = tri[1,i] - biU[1,i-1]*biL[1,i]
    tri[2,i] = tri[2,i] - biU[2,i]*biL[1,i] - biU[1,i]*biL[2,i]
    tri[3,i] = tri[3,i] - biU[2,i+1]*biL[2,i]
  end
  tri[1,n] = tri[1,n] - biU[1,n-1]*biL[1,n]
  tri[2,n] = tri[2,n] - biU[2,n]*biL[1,n] - biU[1,n]*biL[2,n]
end

@inline function mulbiUv!(u,biU,v)
  n = size(biU,2)
  for i=1:n
    u[i] = u[i] + biU[2,i]*v[i] + biU[1,i]*v[i+1]
  end
end

@inline function mulbiLv!(u,biL,v)
  n = size(biL,2)
  u[1] = u[1] + biL[1,1]*v[1]
  for i=2:n
    u[i] = u[i] + biL[2,i-1]*v[i-1] + biL[1,i]*v[i]
  end
  u[n+1] = u[n+1] + biL[2,n]*v[n]
end

@kernel inbounds = true function TriDiagKernel!(tri,@Const(JRhoW),@Const(JWRho),@Const(JWRhoTh),@Const(JRhoThW),fac)

  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  invfac = 1 / fac
  invfac2 = invfac / fac
  if IC <= NumG
    @views @. tri[1,:,IC] = 0
    @views @. tri[2,:,IC] = invfac2
    @views @. tri[3,:,IC] = 0
    @views mulUL!(tri[:,:,IC],JWRho[:,:,IC],JRhoW[:,:,IC])
    @views mulUL!(tri[:,:,IC],JWRhoTh[:,:,IC],JRhoThW[:,:,IC])
  end    
end

@kernel inbounds = true function TriDiagKernel1!(tri,@Const(JRhoW),@Const(JWRho),@Const(JWRhoTh),@Const(JRhoThW),fac)
  Iz,IC, = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]
  Nz = @uniform @ndrange()[1]

  invfac = eltype(tri)(1) / fac
  invfac2 = invfac / fac
  if IC <= NumG
    if Iz == 1  
      tri[2,Iz,IC] = invfac2  - JWRho[2,Iz,IC] * JRhoW[1,Iz,IC] - JWRho[1,Iz,IC] * JRhoW[2,Iz,IC] -
        JWRhoTh[2,1,IC] * JRhoThW[1,Iz,IC] - JWRhoTh[1,Iz,IC] * JRhoThW[2,Iz,IC]
      tri[3,Iz,IC] = - JWRho[2,Iz,IC] * JRhoW[2,Iz,IC] -
        JWRhoTh[2,2,IC] * JRhoThW[2,Iz,IC]
    elseif Iz == Nz  
      tri[1,Iz,IC] = - JWRho[1,Iz-1,IC] * JRhoW[1,Iz,IC] -
        JWRhoTh[1,Iz-1,IC] * JRhoThW[1,Iz,IC]
      tri[2,Iz,IC] = invfac2 - JWRho[2,Iz,IC] * JRhoW[1,Iz,IC] - JWRho[1,Iz,IC] * JRhoW[2,Iz,IC] -
        JWRhoTh[2,Iz,IC] * JRhoThW[1,Iz,IC] - JWRhoTh[1,Iz,IC] * JRhoThW[2,Iz,IC]
    else  
      tri[1,Iz,IC] = - JWRho[1,Iz-1,IC] * JRhoW[1,Iz,IC] -
        JWRhoTh[1,Iz-1,IC] * JRhoThW[1,Iz,IC]
      tri[2,Iz,IC] = invfac2 - JWRho[2,Iz,IC] * JRhoW[1,Iz,IC] - JWRho[1,Iz,IC] * JRhoW[2,Iz,IC] -
       JWRhoTh[2,Iz,IC] * JRhoThW[1,Iz,IC] - JWRhoTh[1,Iz,IC] * JRhoThW[2,Iz,IC]
      tri[3,Iz,IC] = - JWRho[2,Iz+1,IC] * JRhoW[2,Iz,IC] -
        JWRhoTh[2,Iz+1,IC] * JRhoThW[2,Iz,IC]
    end  
  end
end

@kernel inbounds = true function SchurSolveFacKernel!(NumVTr,Nz,k,v,tri,@Const(JRhoW),@Const(JWRho),@Const(JWRhoTh),@Const(JRhoThW),fac)
  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  invfac = 1 / fac
  invfac2 = invfac / fac

  if IC <= NumG
    @views rRho=v[:,IC,1]
    @views rTh=v[:,IC,5]
    @views rw=v[1:Nz-1,IC,4]
    @views sw=k[1:Nz-1,IC,4]
    k[end,IC,4] = 0
    @views @. tri[1,:,IC] = 0
    @views @. tri[2,:,IC] = invfac2
    @views @. tri[3,:,IC] = 0
    @views mulUL!(tri[:,:,IC],JWRho[:,:,IC],JRhoW[:,:,IC])
    @views mulUL!(tri[:,:,IC],JWRhoTh[:,:,IC],JRhoThW[:,:,IC])
    @. rw = invfac * rw
    @views mulbiUv!(rw,JWRho[:,:,IC],rRho)
    @views mulbiUv!(rw,JWRhoTh[:,:,IC],rTh)
    @views triSolve!(sw,tri[:,:,IC],rw)
    @views mulbiLv!(rRho,JRhoW[:,:,IC],sw)
    @views mulbiLv!(rTh,JRhoThW[:,:,IC],sw)
    for iz = 1 : Nz
      k[iz,IC,1] = fac * v[iz,IC,1]
      k[iz,IC,2] = fac * v[iz,IC,2]
      k[iz,IC,3] = fac * v[iz,IC,3]
      k[iz,IC,5] = fac * v[iz,IC,5]
      for iT = 6 : NumVTr
        k[iz,IC,iT] = fac * v[iz,IC,iT]
      end
    end 
  end    
end
@kernel inbounds = true function SchurSolveKernelF!(k,v,@Const(JWRho),@Const(JWRhoTh),fac)
  Iz,IC, = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]
  Nz = @uniform @ndrange()[1]
  if IC <= NumG
    if Iz <= Nz
      v[Iz,IC,4] = v[Iz,IC,4] / fac + JWRho[2,Iz,IC] * v[Iz,IC,1] + JWRho[1,Iz,IC] * v[Iz+1,IC,1] +
        JWRhoTh[2,Iz,IC] * v[Iz,IC,5] + JWRhoTh[1,Iz,IC] * v[Iz+1,IC,5]          
    else
      k[Iz,IC,4] = 0
    end    
  end  
end

@kernel inbounds = true function SchurSolveKernelB!(NumVTr,k,v,@Const(JRhoW),@Const(JRhoThW),fac)
  Iz,IC, = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]
  Nz = @uniform @ndrange()[1]

  if IC <= NumG
    if Iz == 1
      v[Iz,IC,1] += JRhoW[1,Iz,IC] * k[Iz,IC,4] 
      v[Iz,IC,5] += JRhoThW[1,Iz,IC] * k[Iz,IC,4] 
      k[Iz,IC,1] = fac * v[Iz,IC,1]
      k[Iz,IC,2] = fac * v[Iz,IC,2]
      k[Iz,IC,3] = fac * v[Iz,IC,3]
      k[Iz,IC,5] = fac * v[Iz,IC,5]
      for iT = 6 : NumVTr
        k[Iz,IC,iT] = fac * v[Iz,IC,iT]
      end
    elseif Iz == Nz  
      v[Iz,IC,1] += JRhoW[2,Iz,IC] * k[Iz,IC,4] 
      v[Iz,IC,5] += JRhoThW[2,Iz,IC] * k[Iz,IC,4] 
      k[Iz,IC,1] = fac * v[Iz,IC,1]
      k[Iz,IC,2] = fac * v[Iz,IC,2]
      k[Iz,IC,3] = fac * v[Iz,IC,3]
      k[Iz,IC,5] = fac * v[Iz,IC,5]
      for iT = 6 : NumVTr
        k[Iz,IC,iT] = fac * v[Iz,IC,iT]
      end
    else  
      v[Iz,IC,1] += JRhoW[1,Iz,IC] * k[Iz,IC,4] + JRhoW[2,Iz-1,IC] *  k[Iz-1,IC,4] 
      v[Iz,IC,5] += JRhoThW[1,Iz,IC] * k[Iz,IC,4] + JRhoThW[2,Iz-1,IC] *  k[Iz-1,IC,4] 
      k[Iz,IC,1] = fac * v[Iz,IC,1]
      k[Iz,IC,2] = fac * v[Iz,IC,2]
      k[Iz,IC,3] = fac * v[Iz,IC,3]
      k[Iz,IC,5] = fac * v[Iz,IC,5]
      for iT = 6 : NumVTr
        k[Iz,IC,iT] = fac * v[Iz,IC,iT]
      end
    end  
  end
end

@kernel inbounds = true function SchurSolveKernel!(NumVTr,Nz,k,v,tri,@Const(JRhoW),@Const(JWRho),@Const(JWRhoTh),@Const(JRhoThW),fac)
  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    @views triSolve!(k[1:Nz-1,IC,4],tri[:,:,IC],v[1:Nz-1,IC,4])
  end    
end

function SchurSolveGPU!(k,v,J,fac,Cache,Global)
  backend = get_backend(k)
  FT = eltype(k)

  Nz = size(k,1)
  NumG = size(k,2)
  NumVTr = size(k,3) 

  group = (Nz,10)
  ndrange = (Nz,NumG)
# group = (1024)
  groupTri = (64)
  ndrangeTri = (NumG)

  if J.CompTri
    KTriDiagKernel! = TriDiagKernel1!(backend,group)
    KTriDiagKernel!(J.tri,J.JRhoW,J.JWRho,J.JWRhoTh,J.JRhoThW,fac,ndrange=ndrange)
    KernelAbstractions.synchronize(backend)
    J.CompTri = false
  end  
  KSchurSolveKernelF! = SchurSolveKernelF!(backend,group)
  KSchurSolveKernelF!(k,v,J.JWRho,J.JWRhoTh,fac,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  KSchurSolveKernel! = SchurSolveKernel!(backend,groupTri)
  KSchurSolveKernel!(NumVTr,Nz,k,v,J.tri,J.JRhoW,J.JWRho,J.JWRhoTh,J.JRhoThW,fac,ndrange=ndrangeTri)
  KernelAbstractions.synchronize(backend)
  KSchurSolveKernelB! = SchurSolveKernelB!(backend,group)
  KSchurSolveKernelB!(NumVTr,k,v,J.JRhoW,J.JRhoThW,fac,ndrange=ndrange)

end

function SchurSolve!(k,v,J,fac,Cache,Global)
#   sw=(spdiags(repmat(invfac2,n,1),0,n,n)-invfac*JWW-JWRho*JRhoW-JWRhoTh*JRhoThW)\
#     (invfac*rw+JWRho*rRho+JWRhoTh*rTh)
  n1 = size(v,1)
  n2 = size(v,2)
  n = n1 * n2
  tri = J.tri
  JWRho = J.JWRho
  JRhoW = J.JRhoW
  JWRhoTh = J.JWRhoTh
  JRhoThW=J.JRhoThW
  JTrW=J.JTrW
  JWW=J.JWW
  JDiff = J.JDiff
  JAdvC = J.JAdvC
  JAdvF = J.JAdvF
  @views CdTh = Cache.Aux2DG[:,:,1]
  @views CdTr = Cache.Aux2DG[:,:,2:end]
  NumV = Global.Model.NumV
  if size(k,3) > Global.Model.NumV
    NumTr = Global.Model.NumTr
  else
    NumTr = 0
  end  

  invfac=1/fac
  invfac2=invfac/fac

  if Global.Model.JacVerticalDiffusion && Global.Model.JacVerticalAdvection
    for in2=1:n2
      if J.CompTri
        @views @. JAdvC[2,:,in2] += invfac
        @views @. JAdvF[2,:,in2] += invfac
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

      @views rTh=v[1:n1-1,in2,4]
      @views sTh=k[1:n1-1,in2,4]
      @views triSolve!(sTh,JAdvF[:,:,in2],rTh)
      @. rTh = invfac * sTh
    end
    if Global.Model.Equation == "CompressibleMoist"
      for in2=1:n2
        JDiff[2,1,in2] += CdTr[1,in2,1] 
        @views rTh=v[:,in2,NumV+1]
        @views sTh=k[:,in2,NumV+1]
        @views triSolve!(sTh,JDiff[:,:,in2],rTh)
        @. rTh = invfac * sTh
        JDiff[2,1,in2] -= CdTr[1,in2,1] 
      end  
    end    
  elseif Global.Model.JacVerticalAdvection
    for in2=1:n2
      if J.CompTri
        @views @. JAdvC[2,:,in2] += invfac
        @views @. JAdvF[2,:,in2] += invfac
      end
      @views rTh = v[:,in2,1]
      @views sTh = k[:,in2,1]
      @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
      @. rTh = invfac * sTh

      @views rTh = v[:,in2,2]
      @views rTh = v[:,in2,5]
      @views sTh = k[:,in2,5]
      @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
      @. rTh = invfac * sTh

      @views rTh = v[:,in2,2]
      @views sTh = k[:,in2,2]
      @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
      @. rTh = invfac * sTh

      @views rTh=v[:,in2,3]
      @views sTh=k[:,in2,3]
      @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
      @. rTh = invfac * sTh

      @views rTh=v[1:n1-1,in2,4]
      @views sTh=k[1:n1-1,in2,4]
      @views triSolve!(sTh,JAdvF[:,:,in2],rTh)
      @. rTh = invfac * sTh
    end
    if Global.Model.Equation == "CompressibleMoist"
      for in2=1:n2
        @views rTh=v[:,in2,NumV+1]
        @views sTh=k[:,in2,NumV+1]
        @views triSolve!(sTh,JAdvC[:,:,in2],rTh)
        @. rTh = invfac * sTh
      end  
    end    
  end    

  if Global.Model.Equation == "Compressible"
    for in2=1:n2
      @views rRho=v[:,in2,1]
      @views rTh=v[:,in2,5]
      @views rw=v[1:n1-1,in2,4]
      @views sw=k[1:n1-1,in2,4]
      k[n1,in2,4] = 0
      if Global.Model.Damping
        if J.CompTri
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2  - invfac * JWW[1,:,in2]
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoTh[:,:,in2],JRhoThW[:,:,in2])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWRhoTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      else
        if J.CompTri  
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2   
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoTh[:,:,in2],JRhoThW[:,:,in2])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWRhoTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      end
      @views mulbiLv!(rRho,JRhoW[:,:,in2],sw)
      @views mulbiLv!(rTh,JRhoThW[:,:,in2],sw)

      @views @. k[:,in2,1] = fac * rRho
      @views @. k[:,in2,2:3] = fac * v[:,in2,2:3]
      @views @. k[:,in2,5] = fac * rTh
      for iT = 1 : NumTr
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
    for in2=1:n2
      @views rRho=v[:,in2,1]
      @views rTh=v[:,in2,5]
      @views rRhoV=v[:,in2,RhoVPos]
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
          @views mulUL!(tri[:,:,in2],JWRhoTh[:,:,in2],JRhoThW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoV[:,:,in2],JTrW[:,:,in2,RhoVPos])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWRhoV[:,:,in2],rRhoV)
        @views mulbiUv!(rw,JWRhoTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      else
        if J.CompTri  
          @views @. tri[1,:,in2] = 0
          @views @. tri[2,:,in2] = invfac2   
          @views @. tri[3,:,in2] = 0
          @views mulUL!(tri[:,:,in2],JWRho[:,:,in2],JRhoW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoTh[:,:,in2],JRhoThW[:,:,in2])
          @views mulUL!(tri[:,:,in2],JWRhoV[:,:,in2],JTrW[:,:,in2,RhoVPos])
        end
        @. rw = invfac * rw
        @views mulbiUv!(rw,JWRho[:,:,in2],rRho)
        @views mulbiUv!(rw,JWRhoV[:,:,in2],rRhoV)
        @views mulbiUv!(rw,JWRhoTh[:,:,in2],rTh)
        @views triSolve!(sw,tri[:,:,in2],rw)
      end
      @views mulbiLv!(rRho,JRhoW[:,:,in2],sw)
      @views mulbiLv!(rTh,JRhoThW[:,:,in2],sw)

      @views @. k[:,in2,1] = fac * rRho
      @views @. k[:,in2,2:3] = fac * v[:,in2,2:3]
      @views @. k[:,in2,5] = fac * rTh
      for iT = 1 : NumTr
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


