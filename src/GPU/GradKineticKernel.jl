@kernel inbounds = true function GradFullKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),@Const(X),
  @Const(JJ),@Const(M),@Const(Glob),Gravitation)

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  KinF = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+2)
  pCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    pCol[I,J,iz] = p[Iz,ind]  
    RhoCol[I,J,iz+1] = U[Iz,ind,1]
    KinF[I,J,1,iz+1] = eltype(F)(0.5) * (U[Iz,ind,2]^2 + U[Iz,ind,3]^2)
    KinF[I,J,2,iz+1] = KinF[I,J,1,iz+1] + 1/2 * U[Iz,ind,4]^2
    if iz == 1 && Iz == 1
      wCol = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * wCol^2
      KinF[I,J,2,iz] = KinF[I,J,1,iz+1]
    elseif iz == 1
      RhoCol[I,J,1] = U[Iz-1,ind,1]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2
      KinF[I,J,2,iz] =  eltype(F)(0.5) * (U[Iz-1,ind,4]^2 + U[Iz-1,ind,2]^2 + U[Iz-1,ind,3]^2)
    else
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2  
    end
    if iz == ColumnTilesDim && Iz < Nz
      pCol[I,J,iz+1] = p[Iz+1,ind]  
      RhoCol[I,J,iz+2] = U[Iz+1,ind,1]
      KinF[I,J,1,iz+2] = eltype(F)(0.5) * (U[Iz+1,ind,2]^2 + U[Iz+1,ind,3]^2)
      KinF[I,J,1,iz+2] +=  eltype(F)(0.5) * U[Iz,ind,4]^2
    end  
  end

  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz < Nz
    GradZ = eltype(F)(0.5) * (pCol[I,J,iz+1] - pCol[I,J,iz])
    Gradw =  GradZ* (dXdxI[3,3,2,ID,Iz,IF] + dXdxI[3,3,1,ID,Iz+1,IF])
    x = eltype(F)(0.5) * (X[ID,2,1,Iz,IF] + X[ID,2,1,Iz+1,IF])
    y = eltype(F)(0.5) * (X[ID,2,2,Iz,IF] + X[ID,1,2,Iz+1,IF])
    z = eltype(F)(0.5) * (X[ID,2,3,Iz,IF] + X[ID,1,3,Iz+1,IF])
    Grav = Gravitation(x,y,z)
    @atomic :monotonic F[Iz,ind,4] += -(Gradw +
      Grav * (RhoCol[I,J,iz] * JJ[ID,2,Iz,IF] + RhoCol[I,J,iz+1] * JJ[ID,1,Iz+1,IF])) /
      (RhoCol[I,J,iz] * M[Iz,ind,2] + RhoCol[I,J,iz+1] * M[Iz+1,ind,1])
  end

  if Iz <= Nz
    DXpCol = D[I,1] * pCol[1,J,iz]
    DYpCol = D[J,1] * pCol[I,1,iz]
    DXKinF1 = D[I,1] * KinF[1,J,1,iz+1]
    DYKinF1 = D[J,1] * KinF[I,1,1,iz+1]
    DXKinF2 = D[I,1] * KinF[1,J,2,iz+1]
    DYKinF2 = D[J,1] * KinF[I,1,2,iz+1]
    for k = 2 : N
      DXKinF1 += D[I,k] * KinF[k,J,1,iz+1]
      DYKinF1 += D[J,k] * KinF[I,k,1,iz+1]
      DXKinF2 += D[I,k] * KinF[k,J,2,iz+1]
      DYKinF2 += D[J,k] * KinF[I,k,2,iz+1]
      DXpCol += D[I,k] * pCol[k,J,iz]
      DYpCol += D[J,k] * pCol[I,k,iz]
    end
    GraduF1 =
      -(dXdxI[1,1,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,1,1,ID,Iz,IF]  * DYKinF1)
    GradvF1 =
      -(dXdxI[1,2,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,2,1,ID,Iz,IF]  * DYKinF1)
    GraduF2 =
      -(dXdxI[1,1,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,1,2,ID,Iz,IF]  * DYKinF2)
    GradvF2 =
      -(dXdxI[1,2,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,2,2,ID,Iz,IF]  * DYKinF2)
    @views Gradu, Gradv = Grad12(DXpCol,DYpCol,dXdxI[1:2,1:2,:,ID,Iz,IF])

    GradwF1 = eltype(F)(0)
    GradwF2 = eltype(F)(0)
    if Iz > 1 
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+1] - KinF[I,J,2,iz] )  
      GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
      GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
      GradwF1 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,1,ID,Iz,IF]
    end  
    if Iz < Nz
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+2] - KinF[I,J,2,iz+1] )  
      GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
      GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
      GradwF2 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    end  
    GradZ = eltype(F)(0.5) * (KinF[I,J,2,iz+1] - KinF[I,J,1,iz+1])
    GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
    GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
    GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
    GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
    GradwF2 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    GradwF1 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,1,ID,Iz,IF]

    x = eltype(F)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(F)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(F)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    Grav = Gravitation(x,y,z)
    GradZ = -Grav * RhoCol[I,J,iz+1] *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,3,1,ID,Iz,IF] + dXdxI[3,3,2,ID,Iz,IF])
    Gradu += GradZ * (dXdxI[3,1,1,ID,Iz,IF] + dXdxI[3,1,2,ID,Iz,IF])
    Gradv += GradZ * (dXdxI[3,2,1,ID,Iz,IF] + dXdxI[3,2,2,ID,Iz,IF])

    @atomic :monotonic F[Iz,ind,2] += (GraduF1 + GraduF2 - Gradu / RhoCol[I,J,iz+1]) / 
      (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz,ind,3] += (GradvF1 + GradvF2 - Gradv / RhoCol[I,J,iz+1]) / 
      (M[Iz,ind,1] + M[Iz,ind,2])

    if Iz < Nz
      GradZ = eltype(F)(0.5) * (pCol[I,J,iz+1] - pCol[I,J,iz])
      Gradw =  GradZ* (dXdxI[3,3,2,ID,Iz,IF] + dXdxI[3,3,1,ID,Iz+1,IF])
      x = eltype(F)(0.5) * (X[ID,2,1,Iz,IF] + X[ID,2,1,Iz+1,IF])
      y = eltype(F)(0.5) * (X[ID,2,2,Iz,IF] + X[ID,1,2,Iz+1,IF])
      z = eltype(F)(0.5) * (X[ID,2,3,Iz,IF] + X[ID,1,3,Iz+1,IF])
      Grav = Gravitation(x,y,z)
      Gradp = -(Gradw +
        Grav * (RhoCol[I,J,iz+1] * JJ[ID,2,Iz,IF] + RhoCol[I,J,iz+2] * JJ[ID,1,Iz+1,IF])) 
      @atomic :monotonic F[Iz,ind,4] += (GradwF2 + Gradp) / 
        (RhoCol[I,J,iz+2] * M[Iz+1,ind,1] + RhoCol[I,J,iz+1] * M[Iz,ind,2])
    end  
    if Iz > 1  
      @atomic :monotonic F[Iz-1,ind,4] += GradwF1 / 
        (RhoCol[I,J,iz] * M[Iz-1,ind,2] + RhoCol[I,J,iz+1] * M[Iz,ind,1])
    end  
  end

end  

@kernel inbounds = true function RhoGradKinKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  KinF = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+2)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    RhoCol[I,J,iz+1] = U[Iz,ind,1]
    KinF[I,J,1,iz+1] = eltype(F)(0.5) * (U[Iz,ind,2]^2 + U[Iz,ind,3]^2)
    KinF[I,J,2,iz+1] = KinF[I,J,1,iz+1] + 1/2 * U[Iz,ind,4]^2
    if iz == 1 && Iz == 1
      wCol = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * wCol^2
      KinF[I,J,2,iz] = KinF[I,J,1,iz+1]
    elseif iz == 1
      RhoCol[I,J,1] = U[Iz-1,ind,1]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2
      KinF[I,J,2,iz] =  eltype(F)(0.5) * (U[Iz-1,ind,4]^2 + U[Iz-1,ind,2]^2 + U[Iz-1,ind,3]^2)
    else
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2  
    end
    if iz == ColumnTilesDim && Iz < Nz
      RhoCol[I,J,iz+2] = U[Iz+1,ind,1]
      KinF[I,J,1,iz+2] = eltype(F)(0.5) * (U[Iz+1,ind,2]^2 + U[Iz+1,ind,3]^2)
      KinF[I,J,1,iz+2] +=  eltype(F)(0.5) * U[Iz,ind,4]^2
    end  
  end

  @synchronize

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    DXKinF1 = D[I,1] * KinF[1,J,1,iz+1]
    DYKinF1 = D[J,1] * KinF[I,1,1,iz+1]
    DXKinF2 = D[I,1] * KinF[1,J,2,iz+1]
    DYKinF2 = D[J,1] * KinF[I,1,2,iz+1]
    for k = 2 : N
      DXKinF1 += D[I,k] * KinF[k,J,1,iz+1]
      DYKinF1 += D[J,k] * KinF[I,k,1,iz+1]
      DXKinF2 += D[I,k] * KinF[k,J,2,iz+1]
      DYKinF2 += D[J,k] * KinF[I,k,2,iz+1]
    end
    GraduF1 =
      -(dXdxI[1,1,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,1,1,ID,Iz,IF]  * DYKinF1)
    GradvF1 =
      -(dXdxI[1,2,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,2,1,ID,Iz,IF]  * DYKinF1)
    GraduF2 =
      -(dXdxI[1,1,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,1,2,ID,Iz,IF]  * DYKinF2)
    GradvF2 =
      -(dXdxI[1,2,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,2,2,ID,Iz,IF]  * DYKinF2)

    GradwF1 = eltype(F)(0)
    GradwF2 = eltype(F)(0)
    if Iz > 1 
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+1] - KinF[I,J,2,iz] )  
      GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
      GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
      GradwF1 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,1,ID,Iz,IF]
    end  
    if Iz < Nz
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+2] - KinF[I,J,2,iz+1] )  
      GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
      GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
      GradwF2 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    end  
    GradZ = eltype(F)(0.5) * (KinF[I,J,2,iz+1] - KinF[I,J,1,iz+1])
    GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
    GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
    GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
    GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
    GradwF2 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    GradwF1 += -RhoCol[I,J,iz+1] * GradZ * dXdxI[3,3,1,ID,Iz,IF]

    @atomic :monotonic F[Iz,ind,2] += (GraduF1 + GraduF2) / (M[Iz,ind,1] + M[Iz,ind,2])
    @atomic :monotonic F[Iz,ind,3] += (GradvF1 + GradvF2) / (M[Iz,ind,1] + M[Iz,ind,2])

    if Iz < Nz
      @atomic :monotonic F[Iz,ind,4] += GradwF2 / 
        (RhoCol[I,J,iz+2] * M[Iz+1,ind,1] + RhoCol[I,J,iz+1] * M[Iz,ind,2])
    end  
    if Iz > 1  
      @atomic :monotonic F[Iz-1,ind,4] += GradwF1 / 
        (RhoCol[I,J,iz] * M[Iz-1,ind,2] + RhoCol[I,J,iz+1] * M[Iz,ind,1])
    end  
  end

end  

@kernel inbounds = true function RhoGradKinEDMFKernel!(F,@Const(U),@Const(w),@Const(Rho),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF,IE = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]
  NE = @uniform @ndrange()[5]

  KinF = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    RhoCol[I,J,iz] = Rho[Iz,ind,IE]
    KinF[I,J,1,iz+1] = eltype(F)(0.5) * (U[Iz,ind,2]^2 + U[Iz,ind,3]^2)
    KinF[I,J,2,iz+1] = KinF[I,J,1,iz+1] + 1/2 * w[Iz,ind,IE]^2
    if iz == 1 && Iz == 1
      wCol = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * wCol^2
      KinF[I,J,2,iz] = KinF[I,J,1,iz+1]
    elseif iz == 1
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * w[Iz-1,ind,IE]^2
      KinF[I,J,2,iz] =  eltype(F)(0.5) * (w[Iz-1,ind,IE]^2 + U[Iz-1,ind,2]^2 + U[Iz-1,ind,3]^2)
    else
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * w[Iz-1,ind,IE]^2  
    end
    if iz == ColumnTilesDim && Iz < Nz
      RhoCol[I,J,iz+1] = Rho[Iz+1,ind,IE]
      KinF[I,J,1,iz+2] = eltype(F)(0.5) * (U[Iz+1,ind,2]^2 + U[Iz+1,ind,3]^2)
      KinF[I,J,1,iz+2] +=  eltype(F)(0.5) * w[Iz,ind,IE]^2
    end  
  end

  @synchronize


  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    GradwF1 = eltype(F)(0)
    GradwF2 = eltype(F)(0)
    RhoCol = Rho[Iz,ind,IE]
    if Iz > 1 
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+1] - KinF[I,J,2,iz] )  
      GradwF1 += -RhoCol * GradZ * dXdxI[3,3,1,ID,Iz,IF]
    end  
    if Iz < Nz
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+2] - KinF[I,J,2,iz+1] )  
      GradwF2 += -RhoCol * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    end  
    GradZ = eltype(F)(0.5) * (KinF[I,J,2,iz+1] - KinF[I,J,1,iz+1])
    GradwF2 += -RhoCol * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    GradwF1 += -RhoCol * GradZ * dXdxI[3,3,1,ID,Iz,IF]


    if Iz < Nz
      @atomic :monotonic F[Iz,ind,4] += GradwF2 / 
      (RhoCol[I,J,iz] * M[Iz,ind,2] + RhoCol[I,J,iz+1] * M[Iz+1,ind,1])
    end  
    if Iz > 1  
      @atomic :monotonic F[Iz-1,ind,4] += GradwF1 / 
        (RhoCol[I,J,iz] * M[Iz,ind,1] + RhoCol[I,J,iz-1] * M[Iz-1,ind,2])
    end  
  end

end  
