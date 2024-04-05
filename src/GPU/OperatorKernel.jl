#=
@kernel inbounds = true function GradKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(MRho),@Const(Glob),Phys)

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  Pres = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    Pres[I,J,1,iz+1] = p[Iz,ind]
    Pres[I,J,2,iz+1] = p[Iz,ind]
  end
  if iz == ColumnTilesDim && Iz < Nz
    Pres[I,J,iz+1] = p[Iz+1,ind]
  end  

  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    DXPres = D[I,1] * Pres[1,J,iz]
    DYPres = D[J,1] * Pres[I,1,iz]
    for k = 2 : N
      DXPres += D[I,k] * Pres[k,J,iz]
      DYPres += D[J,k] * Pres[I,k,iz]
    end
    @views Gradu, Gradv = Grad12(DXPres,DYPres,dXdxI[1:2,1:2,:,ID,Iz,IF]) 

    GradZ = -Phys.Grav * U[Iz,ind,1] *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,3,1,ID,Iz,IF] + dXdxI[3,3,2,ID,Iz,IF])
    Gradu += GradZ * (dXdxI[3,1,1,ID,Iz,IF] + dXdxI[3,1,2,ID,Iz,IF])
    Gradv += GradZ * (dXdxI[3,2,1,ID,Iz,IF] + dXdxI[3,2,2,ID,Iz,IF])
    @atomic :monotonic F[Iz,ind,2] += -Gradu / M[Iz,ind] / U[Iz,ind,1]
    @atomic :monotonic F[Iz,ind,3] += -Gradv / M[Iz,ind] / U[Iz,ind,1]
  end  

  if Iz < Nz
    GradZ = eltype(F)(0.5) * (Pres[I,J,iz+1] - Pres[I,J,iz])  
    Gradw =  GradZ* (dXdxI[3,3,2,ID,Iz,IF] + dXdxI[3,3,1,ID,Iz+1,IF])
    @atomic :monotonic F[Iz,ind,4] += -(Gradw +
      Phys.Grav * (U[Iz,ind,1] * JJ[ID,2,Iz,IF] + U[Iz+1,ind,1] * JJ[ID,1,Iz+1,IF])) /
      MRho[Iz,ind]
  end      
   
end

@kernel inbounds = true function GradDeepKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(MRho),@Const(Glob),Phys)

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  Pres = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    Pres[I,J,iz] = p[Iz,ind]
  end
  if iz == ColumnTilesDim && Iz < Nz
    Pres[I,J,iz+1] = p[Iz+1,ind]
  end  

  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    DXPres = D[I,1] * Pres[1,J,iz]
    DYPres = D[J,1] * Pres[I,1,iz]
    for k = 2 : N
      DXPres += D[I,k] * Pres[k,J,iz]
      DYPres += D[J,k] * Pres[I,k,iz]
    end
    @views Gradu, Gradv = Grad12(DXPres,DYPres,dXdxI[1:2,1:2,:,ID,Iz,IF]) 

    GradZ = -Phys.Grav * U[Iz,ind,1] *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,3,1,ID,Iz,IF] + dXdxI[3,3,2,ID,Iz,IF])
    Gradu += GradZ * (dXdxI[3,1,1,ID,Iz,IF] + dXdxI[3,1,2,ID,Iz,IF])
    Gradv += GradZ * (dXdxI[3,2,1,ID,Iz,IF] + dXdxI[3,2,2,ID,Iz,IF])
    @atomic :monotonic F[Iz,ind,2] += -Gradu / M[Iz,ind] / U[Iz,ind,1]
    @atomic :monotonic F[Iz,ind,3] += -Gradv / M[Iz,ind] / U[Iz,ind,1]
  end  

  if Iz < Nz
    GradZ = eltype(F)(0.5) * (Pres[I,J,iz+1] - Pres[I,J,iz])  
    Gradw =  GradZ* (dXdxI[3,3,2,ID,Iz,IF] + dXdxI[3,3,1,ID,Iz+1,IF])
    x = X[ID,2,1,Iz,IF] 
    y = X[ID,2,2,Iz,IF]
    z = X[ID,2,3,Iz,IF]
    r = sqrt(x^2 + y^2 + z^2)
    @atomic :monotonic F[Iz,ind,4] += -(Gradw + 
      Phys.Grav * (Phys.RadEarth / r)^2 * (U[Iz,ind,1] * JJ[ID,2,Iz,IF] + U[Iz+1,ind,1] * JJ[ID,1,Iz+1,IF])) /
      MRho[Iz,ind]
  end      
   
end
=#

@kernel inbounds = true function DivRhoGradKernel!(F,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim)
  FCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  if Iz <= Nz
    cCol[I,J,iz] = U[Iz,ind,5] / U[Iz,ind,1]
    FCol[I,J,iz] = 0.0
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    Dxc = D[I,1] * cCol[1,J,iz]
    Dyc = D[J,1] * cCol[I,1,iz]
    for k = 2 : N
      Dxc = Dxc + D[I,k] * cCol[k,J,iz]
      Dyc = Dyc + D[J,k] * cCol[I,k,iz] 
    end
    @views (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @atomic :monotonic FCol[k,J,iz] += DW[k,I] * tempx
      @atomic :monotonic FCol[I,k,iz] += DW[k,J] * tempy
    end
  end

  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz
    @atomic :monotonic F[Iz,ind,5] += FCol[I,J,iz] / M[Iz,ind]
  end
end


@kernel inbounds = true function DivRhoTrCentralKernel!(F,@Const(c),@Const(uC),@Const(vC),@Const(w),
  @Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  uCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  vCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  wCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+1)
  FCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  if Iz <= Nz
    wCol[I,J,iz+1] = w[Iz,ind]
    cCol[I,J,iz+1] = c[Iz,ind]
    uCol[I,J,iz+1] = uC[Iz,ind]
    vCol[I,J,iz+1] = vC[Iz,ind]
    FCol[I,J,iz+1] = 0
    if iz == 1 && Iz > 1
      cCol[I,J,1] = c[Iz-1,ind]
      uCol[I,J,1] = uC[Iz-1,ind]
      vCol[I,J,1] = vC[Iz-1,ind]
      wCol[I,J,1] = w[Iz,ind]
      FCol[I,J,1] = 0
    elseif iz == 1 && Iz == 1
      cCol[I,J,1] = c[1,ind]
      wCol[I,J,1] = 0
      FCol[I,J,1] = 0
    end
    if iz == ColumnTilesDim && Iz < Nz
      cCol[I,J,ColumnTilesDim+2] = c[Iz+1,ind]
      uCol[I,J,ColumnTilesDim+2] = uC[Iz+1,ind]
      vCol[I,J,ColumnTilesDim+2] = vC[Iz+1,ind]
      FCol[I,J,ColumnTilesDim+2] = 0
    elseif iz == ColumnTilesDim && Iz == Nz
      cCol[I,J,ColumnTilesDim+2] = c[Nz,ind]
      FCol[I,J,ColumnTilesDim+2] = 0
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    wCon = dXdxI[I,J,2,Iz,3,1,IF] * uCol[I,J,iz+1] + 
      dXdxI[I,J,2,Iz,3,2,IF] * vCol[I,J,iz+1] + 
      dXdxI[I,J,1,Iz+1,3,1,IF] * uCol[I,J,iz+2] + 
      dXdxI[I,J,2,Iz+1,3,2,IF] * vCol[I,J,iz+2] + 
      (dXdxI[I,J,2,Iz,3,3,IF] + dXdxI[I,J,1,Iz+1,3,3,IF]) * wCol[I,J,iz+1]
    cF = (JJ[ID,2,Iz,IF] * cCol[I,J,iz+1] + JJ[ID,1,Iz+1,IF] * cCol[I,J,iz+2]) /
      (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    Flux = eltype(F)(0.5) * wCon * cF
    @atomic :monotonic FCol[I,J,iz+1] += -Flux
    @atomic :monotonic FCol[I,J,iz+2] += Flux
  end 

  if Iz <= Nz
    tempx = -cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,1,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,1,3,IF] * wCol[I,J,iz+1])
    tempy = -cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,2,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,2,3,IF] * wCol[I,J,iz+1])
    for k = 1 : N
      @atomic :monotonic FCol[k,J,iz+1] += D[k,I] * tempx
      @atomic :monotonic FCol[I,k,iz+1] += D[k,J] * tempy
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz 
    @atomic :monotonic F[Iz,ind] += FCol[I,J,iz+1] / M[Iz,ind]
    if iz == 1 && Iz >  1
      @atomic :monotonic F[Iz-1,ind] += FCol[I,J,iz] / M[Iz-1,ind]
    end
    if iz == ColumnTilesDim && Iz <  Nz
      @atomic :monotonic F[Iz+1,ind] += FCol[I,J,iz+2] / M[Iz+1,ind]
    end
  end
end

@kernel inbounds = true function DivRhoTrUpwindKernel!(F,@Const(c),@Const(Rho),@Const(uC),@Const(vC),@Const(w),
  @Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim+1)
  FCol = @localmem eltype(F) (N,N, ColumnTilesDim+2)
  if Iz <= Nz
    wCol[I,J,iz+1] = w[Iz,ind]
    RhoCol[I,J,iz+1] = Rho[Iz,ind]
    cCol[I,J,iz+1] = c[Iz,ind] / RhoCol[I,J,iz+1]
    uCol[I,J,iz+1] = uC[Iz,ind]
    vCol[I,J,iz+1] = vC[Iz,ind]
    FCol[I,J,iz+1] = 0
    if iz == 1 && Iz > 1
      RhoCol[I,J,1] = Rho[Iz-1,ind]
      cCol[I,J,1] = c[Iz-1,ind] / RhoCol[I,J,1]
      uCol[I,J,1] = uC[Iz-1,ind]
      vCol[I,J,1] = vC[Iz-1,ind]
      wCol[I,J,1] = w[Iz,ind]
      FCol[I,J,1] = 0
    elseif iz == 1 && Iz == 1
      RhoCol[I,J,1] = Rho[1,ind]
      cCol[I,J,1] = c[1,ind] / RhoCol[I,J,1]
      wCol[I,J,1] = 0
      FCol[I,J,1] = 0
    end
    if iz == ColumnTilesDim && Iz < Nz
      RhoCol[I,J,ColumnTilesDim+2] = Rho[Iz+1,ind]
      cCol[I,J,ColumnTilesDim+2] = c[Iz+1,ind] / RhoCol[I,J,ColumnTilesDim+2]
      uCol[I,J,ColumnTilesDim+2] = uC[Iz+1,ind]
      vCol[I,J,ColumnTilesDim+2] = vC[Iz+1,ind]
      FCol[I,J,ColumnTilesDim+2] = 0
    elseif iz == ColumnTilesDim && Iz == Nz
      RhoCol[I,J,ColumnTilesDim+2] = Rho[Nz,ind]
      cCol[I,J,ColumnTilesDim+2] = c[Nz,ind] / RhoCol[I,J,ColumnTilesDim+2]
      FCol[I,J,ColumnTilesDim+2] = 0
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    wCon = RhoCol[I,J,iz+1] * (dXdxI[I,J,2,Iz,3,1,IF] * uCol[I,J,iz+1] + 
      dXdxI[I,J,2,Iz,3,2,IF] * vCol[I,J,iz+1] + dXdxI[I,J,2,Iz,3,3,IF] * wCol[I,J,iz+1]) +
      RhoCol[I,J,iz+2] * (dXdxI[I,J,1,Iz+1,3,1,IF] * uCol[I,J,iz+2] + 
      dXdxI[I,J,2,Iz+1,3,2,IF] * vCol[I,J,iz+2] + dXdxI[I,J,1,Iz+1,3,3,IF] * wCol[I,J,iz+1])
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    Flux = 0.25 * ((abs(wCon) + wCon) * cL + (-abs(wCon) + wCon) * cR)
    @atomic :monotonic FCol[I,J,iz+1] += -Flux
    @atomic :monotonic FCol[I,J,iz+2] += Flux
  end 

  if Iz <= Nz
    tempx = -RhoCol[I,J,iz+1] * cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,1,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,1,3,IF] * wCol[I,J,iz+1])
    tempy = -RhoCol[I,J,iz+1] * cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,2,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,2,3,IF] * wCol[I,J,iz+1])
    for k = 1 : N
      @atomic :monotonic FCol[k,J,iz+1] += D[k,I] * tempx
      @atomic :monotonic FCol[I,k,iz+1] += D[k,J] * tempy
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]
  if Iz <= Nz 
    @atomic :monotonic F[Iz,ind] += FCol[I,J,iz+1] / M[Iz,ind]
    if iz == 1 && Iz >  1
      @atomic :monotonic F[Iz-1,ind] += FCol[I,J,iz] / M[Iz-1,ind]
    end
    if iz == ColumnTilesDim && Iz <  Nz
      @atomic :monotonic F[Iz+1,ind] += FCol[I,J,iz+2] / M[Iz+1,ind]
    end
  end
end

@kernel inbounds = true function DivRhoThUpwind3Kernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  if Iz <= Nz
    cCol[I,J,iz+1] = U[Iz,ind,5] / U[Iz,ind,1]
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = U[Izm1,ind,5] / U[Izm1,ind,1]
#   if Iz == 1
#     cCol[I,J,iz] = eltype(F)(2) * cCol[I,J,iz] - U[2,ind,5] / U[2,ind,1]  
#   end  
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = U[Izp1,ind,5] / U[Izp1,ind,1]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = U[Izp2,ind,5] / U[Izp2,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3W(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF],JJ[ID,:,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(F)(0.5) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic F[Iz,ind,5] += -Flux / M[Iz,ind]
    @atomic :monotonic F[Iz+1,ind,5] += Flux / M[Iz+1,ind]
    Flux = wCon
    @atomic :monotonic F[Iz,ind,1] += -Flux / M[Iz,ind]
    @atomic :monotonic F[Iz+1,ind,1] += Flux / M[Iz+1,ind]
  end 

  if Iz <= Nz
    DivRho = D[I,1] * uConCol[1,J,iz] 
    DivRho += D[J,1] * vConCol[I,1,iz] 
    DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] 
    DivRhoTr += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRho += D[I,k] * uConCol[k,J,iz] 
      DivRho += D[J,k] * vConCol[I,k,iz] 
      DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] 
      DivRhoTr += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic F[Iz,ind,1] += DivRho / M[Iz,ind]
    @atomic :monotonic F[Iz,ind,5] += DivRhoTr / M[Iz,ind]
  end
end

@kernel inbounds = true function DivRhoKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  if Iz <= Nz
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Flux = eltype(F)(0.5) * wCon
    @atomic :monotonic F[Iz,ind,1] += -Flux / M[Iz,ind]
    @atomic :monotonic F[Iz+1,ind,1] += Flux / M[Iz+1,ind]
  end 

  if Iz <= Nz
    DivRho = D[I,1] * uConCol[1,J,iz] 
    DivRho += D[J,1] * vConCol[I,1,iz] 
    for k = 2 : N
      DivRho += D[I,k] * uConCol[k,J,iz] 
      DivRho += D[J,k] * vConCol[I,k,iz] 
    end
    @atomic :monotonic F[Iz,ind,1] += DivRho / M[Iz,ind]
  end
end

@kernel inbounds = true function DivRhoTrUpwind3Kernel!(FTr,@Const(Tr),@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  uConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind] / U[Iz,ind,1]
    @views (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    uConCol[I,J,iz] = uCon
    vConCol[I,J,iz] = vCon
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind] / U[Izm1,ind,1]
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind] / U[Izp1,ind,1]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind] / U[Izp2,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind] += -Flux / M[Iz,ind]
    @atomic :monotonic FTr[Iz+1,ind] += Flux / M[Iz+1,ind]
  end 

  if Iz <= Nz
    DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz+1] 
    DivRhoTr += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz+1]
    for k = 2 : N
      DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz+1] 
      DivRhoTr += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz+1]
    end
    @atomic :monotonic FTr[Iz,ind] += DivRhoTr / M[Iz,ind]
  end
end

@kernel inbounds = true function AdvectionTrUpwind3Kernel!(FTr,@Const(Tr),@Const(U),@Const(w),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(FTr) (N,N, ColumnTilesDim+3)
  if Iz <= Nz
    cCol[I,J,iz+1] = Tr[Iz,ind]
  end
  if iz == 1
    Izm1 = max(Iz - 1,1)
    cCol[I,J,iz] = Tr[Izm1,ind] 
  end
  if iz == ColumnTilesDim || Iz == Nz
    Izp1 = min(Iz + 1,Nz)
    cCol[I,J,iz+2] = Tr[Izp1,ind]
    Izp2 = min(Iz + 2,Nz)
    cCol[I,J,iz+3] = Tr[Izp2,ind]
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cLL = cCol[I,J,iz]
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    cRR = cCol[I,J,iz+3]

    @views wCon = Contra3(U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      w[Iz,ind],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = eltype(FTr)(0.25) * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic FTr[Iz,ind] += (-Flux + wCon * cL) / M[Iz,ind] 
    @atomic :monotonic FTr[Iz+1,ind] += (Flux - wCon * cR) / M[Iz+1,ind]
  end 

  if Iz <= Nz
    GradxTr = D[I,1] * cCol[1,J,iz+1] 
    GradyTr = D[J,1] * cCol[I,1,iz+1]
    for k = 2 : N
      GradxTr += D[I,k] * cCol[k,J,iz+1] 
      GradyTr += D[J,k] * cCol[I,k,iz+1]
    end
    @views (uCon, vCon) = Contra12(U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @atomic :monotonic FTr[Iz,ind] += -(GradxTr * uCon + GradyTr * vCon) / M[Iz,ind]
  end
end

#=
@kernel inbounds = true function DivRhoTrUpwind3Kernel!(F,@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),Koeff)

# Combines hyperViscosity with advection for tracers

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  CacheCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FTrCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FRhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  if Iz <= Nz
    CacheCol[I,J,iz] = Cache[Iz,ind]
    wCol[I,J,iz] = U[Iz,ind,4]
    RhoCol[I,J,iz] = U[Iz,ind,1]
    cCol[I,J,iz] = U[Iz,ind,5] / RhoCol[I,J,iz]
    uCol[I,J,iz] = U[Iz,ind,2]
    vCol[I,J,iz] = U[Iz,ind,3]
    FRhoCol[I,J,iz] = 0
    FTrCol[I,J,iz] = 0
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz < Nz 
    cL = cCol[I,J,iz]
    cR = cCol[I,J,iz+1]
    if iz > 1
      cLL = cCol[I,J,iz-1]
    else
      Izm1 = max(Iz - 1,1)
      cLL = U[Izm1,ind,5] / U[Izm1,ind,1]
    end
    if iz < ColumnTilesDim - 1
      cRR = cCol[I,J,iz+2]
    else
      Izp2 = min(Iz + 2, Nz)
      cRR = U[Izp2,ind,5] / U[Izp2,ind,1]
    end

    @views wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = 0.25 * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @atomic :monotonic F[Iz,ind,5] += -Flux / M[Iz,ind]
    @atomic :monotonic F[Iz+1,ind,5] += Flux / M[Iz+1,ind]
    Flux = eltype(F)(0.5) * wCon
    @atomic :monotonic F[Iz,ind,1] += -Flux / M[Iz,ind]
    @atomic :monotonic F[Iz+1,ind,1] += Flux / M[Iz+1,ind]
  end 

  if Iz <= Nz
    Dxc = D[I,1] * CacheCol[1,J,iz]
    Dyc = D[J,1] * CacheCol[I,1,iz]
    for k = 2 : N
      Dxc += D[I,k] * CacheCol[k,J,iz]
      Dyc += D[J,k] * CacheCol[I,k,iz]
    end
    
    @views (GradDx, GradDy) = Grad12(RhoCol[I,J,iz],Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views (tempx, tempy) = Contra12(-Koeff,GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @atomic :monotonic FTrCol[k,J,iz] += DW[k,I] * tempx
      @atomic :monotonic FTrCol[I,k,iz] += DW[k,J] * tempy
    end

    @views (tempxRho, tempyRho) = Contra12(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @atomic :monotonic FRhoCol[k,J,iz] += D[k,I] * tempxRho
      @atomic :monotonic FRhoCol[I,k,iz] += D[k,J] * tempyRho
    end
    tempxTr = tempxRho * cCol[I,J,iz]
    tempyTr = tempyRho * cCol[I,J,iz]
    for k = 1 : N
      @atomic :monotonic FTrCol[k,J,iz] += D[k,I] * tempxTr
      @atomic :monotonic FTrCol[I,k,iz] += D[k,J] * tempyTr
    end
  end
  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz 
    @atomic :monotonic F[Iz,ind,1] += FRhoCol[I,J,iz] / M[Iz,ind]
    @atomic :monotonic F[Iz,ind,5] += FTrCol[I,J,iz] / M[Iz,ind]
  end
end
=#


@inline function Contra12(Rho,u,v,dXdxI)
  uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * v)
  vCon = Rho * ((dXdxI[2,1,1] + dXdxI[2,1,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v)
  return uCon, vCon
end

@inline function Contra12(u,v,dXdxI)
  uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * v
  vCon = (dXdxI[2,1,1] + dXdxI[2,1,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  return uCon, vCon
end

@inline function Grad12(Rho,u,v,dXdxI,J)
  uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v) / (J[1] + J[2])
  vCon = Rho * ((dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v) / (J[1] + J[2])
  return uCon, vCon
end

@inline function Grad12(u,v,dXdxI,J)
  uCon = ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v) / (J[1] + J[2])
  vCon = ((dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v) / (J[1] + J[2])
  return uCon, vCon
end

@inline function Grad12(u,v,dXdxI)
  uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v
  vCon = (dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  return uCon, vCon
end

@inline function Grad3(u,v,dXdxI)
  wCon1 = dXdxI[1,3,1] * u + dXdxI[2,3,1] * v
  wCon2 = dXdxI[1,3,2] * u + dXdxI[2,3,2] * v
  return wCon1, wCon2
end

@inline function Curl12(u,v,dXdxI)
  uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * v -
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * u
  vCon = (dXdxI[2,1,1] + dXdxI[2,1,2]) * v -
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * u
  return uCon, vCon
end

@inline function Rot12(u,v,dXdxI)
  uCon = (dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  vCon = -(dXdxI[1,1,1] + dXdxI[1,1,2]) * u -
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v
  return uCon, vCon
end

@inline function Contra3(Rho,u,v,w,dXdxI)
  wCon = Rho[1] * (dXdxI[1,2,1] * u[1] + dXdxI[2,2,1] * v[1] + dXdxI[3,2,1] * w) + 
    Rho[2] * (dXdxI[1,1,2] * u[2] + dXdxI[2,1,2] * v[2] + dXdxI[3,1,2] * w)
end

@inline function Contra3W(Rho,u,v,w,dXdxI,J)
  wCon = eltype(Rho)(2) * (J[2,1] * Rho[1] * (dXdxI[1,2,1] * u[1] + dXdxI[2,2,1] * v[1] + dXdxI[3,2,1] * w) + 
    J[1,2] * Rho[2] * (dXdxI[1,1,2] * u[2] + dXdxI[2,1,2] * v[2] + dXdxI[3,1,2] * w)) /
    (J[2,1] + J[1,2])
end
  
@inline function RecU41(cLL,cL,cR,cRR,JLL,JL,JR,JRR)

  kR = (JL / (JL + JR)) * ((JLL + JL) / (JLL + JL + JR))
  kL = -(JL / (JLL + JL)) * (JR / (JLL + JL + JR))
  cFL = kL * cLL + (1 - kL - kR)*cL + kR * cR

  kL = (JR / (JR + JL)) * ((JRR + JR)/(JL + JR + JRR))
  kR = -(JR /(JRR + JR)) *(JL /(JL + JR + JRR))
  cFR = kL * cL + (1 - kL - kR) * cR + kR * cRR
 
  return (cFL,cFR)
end

@inline function RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR)

  kR = (JL / (JL + JR)) * ((JLL + JL) / (JLL + JL + JR))
  kL = (JL / (JLL + JL)) * (JR / (JLL + JL + JR))
  #r = (cL - cLL) / (cR - cL + sign(cR- cL)*1.e-20 + 1.e-20)
  r = (cL - cLL) / (cR - cL + 1.e-20)
  cFL = cL + max(eltype(cLL)(0), min(r, min(kL * r + kR, eltype(cLL)(1)))) * (cR - cL) 

  kL = (JR / (JR + JL)) * ((JRR + JR)/(JL + JR + JRR))
  kR = (JR /(JRR + JR)) *(JL /(JL + JR + JRR))
# r = (cR - cRR) / (cL - cR  + sign(cL- cR)*1.e-20 + 1.e-20)
  r = (cR - cRR) / (cL - cR  + 1.e-20)
  cFR = cR + max(eltype(cLL)(0), min(r, min(kR * r + kL, eltype(cLL)(1)))) * (cL - cR) 
 
  return (cFL,cFR)
end

@kernel inbounds = true function DampKernel!(Damp,F,U,zP)
  Iz,IC = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    Fu,Fv,Fw = Damp(zP[Iz,IC],view(U,Iz,IC,1:5))
    F[Iz,IC,2] += Fu
    F[Iz,IC,3] += Fv
    F[Iz,IC,4] += Fw
  end
end  

@kernel inbounds = true function ForceKernel!(Force,F,U,p,xS)
  Iz,IC = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    FRho,Fu,Fv,Fw,FRhoTh = Force(view(U,Iz,IC,1:5),p[Iz,IC],xS[2,IC])
    F[Iz,IC,1] += FRho
    F[Iz,IC,2] += Fu
    F[Iz,IC,3] += Fv
    F[Iz,IC,4] += Fw
    F[Iz,IC,5] += FRhoTh
  end
end  

@kernel inbounds = true function MicrophysicsKernel!(Source,F,U,p)
  Iz,IC = @index(Global, NTuple)
  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    FRho,FRhoTh,FRhoV,FRhoC = Source(view(U,Iz,IC,:),p[Iz,IC])
    F[Iz,IC,1] += FRho
    F[Iz,IC,5] += FRhoTh
    F[Iz,IC,6] += FRhoV
    F[Iz,IC,7] += FRhoC
  end
end  

@kernel inbounds = true function VerticalDiffusionScalarKernel!(FTr,@Const(Tr),@Const(Rho),@Const(K),
  @Const(dz))
  iz,iC = @index(Local, NTuple)
  Iz,IC = @index(Global, NTuple)

  Nz = @uniform @ndrange()[1]
  NumG = @uniform @ndrange()[2]

  if Iz < Nz && IC <= NumG
    grad = eltype(FTr)(2) * K[Iz,IC] * (Tr[Iz+1,IC] / Rho[Iz+1,IC] - 
      Tr[Iz,IC] / Rho[Iz,IC]) / (dz[Iz+1,IC] + dz[Iz,IC])
    @atomic :monotonic FTr[Iz,IC] +=  grad / dz[Iz,IC]
    @atomic :monotonic FTr[Iz+1,IC] +=  -grad / dz[Iz+1,IC]
  end  
end  

@kernel inbounds = true function SurfaceFluxScalarsKernel!(F,@Const(U),@Const(p),@Const(TSurf),@Const(RhoVSurf),@Const(uStar),
  @Const(CT),@Const(CH),@Const(dz),@Const(Phys))

  IC, = @index(Global, NTuple)

  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    Rho = U[1,IC,1]
    RhoTh = U[1,IC,5]
    RhoV = U[1,IC,6]
    dz1 = dz[1,IC]
    RhoD = Rho - RhoV
    Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV
    T = p[1,IC] / Rm
    LatFlux = -CT[IC] * uStar[IC] * (RhoV - RhoVSurf[IC]) 
    SensFlux = -CH[IC] * uStar[IC] * (T - TSurf[IC]) 
    FRho = LatFlux
    FRhoV = LatFlux
    PrePi=(p[1,IC] / Phys.p0)^(Rm / Cpml)
    FRhoTh = RhoTh * (SensFlux / T + ((Phys.Rv / Rm) - eltype(F)(1) / Rho - 
      log(PrePi)*(Phys.Rv / Rm - Phys.Cpv / Cpml)) *  LatFlux)
    F[1,IC,1] += FRho / dz1 
    F[1,IC,5] += FRhoTh  / dz1
    F[1,IC,6] += FRhoV  / dz1
#   @views SurfaceFluxScalars(FU[1,IC,:],U[1,IC,:],p[1,IC],CT[IC],CH[IC],uStar[iC])
  end  
end

