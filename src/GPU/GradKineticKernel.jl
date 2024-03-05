@kernel inbounds = true function RhoGradKinKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(MRho),@Const(Glob))

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  KinF = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)

  ID = I + (J - 1) * N
  ind = Glob[ID,IF]

  if Iz <= Nz
    KinF[I,J,1,iz+1] = eltype(F)(0.5) * (U[Iz,ind,2]^2 + U[Iz,ind,3]^2)
    KinF[I,J,2,iz+1] = KinF[I,J,1,iz+1] + 1/2 * U[Iz,ind,4]^2
    if iz == 1 && Iz == 1
      wCol = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] +
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * wCol^2
      KinF[I,J,2,iz] = KinF[I,J,1,iz+1]
    elseif iz == 1
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2
      KinF[I,J,2,iz] =  eltype(F)(0.5) * (U[Iz-1,ind,4]^2 + U[Iz-1,ind,2]^2 + U[Iz-1,ind,3]^2)
    else
      KinF[I,J,1,iz+1] +=  eltype(F)(0.5) * U[Iz-1,ind,4]^2  
    end
    if iz == ColumnTilesDim && Iz < Nz
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
    RhoCol = U[Iz,ind,1]
    if Iz > 1 
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+1] - KinF[I,J,2,iz] )  
      GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
      GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
      GradwF1 += -RhoCol * GradZ * dXdxI[3,3,1,ID,Iz,IF]
    end  
    if Iz < Nz
      GradZ = eltype(F)(0.5) * (KinF[I,J,1,iz+2] - KinF[I,J,2,iz+1] )  
      GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
      GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
      GradwF2 += -RhoCol * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    end  
    GradZ = eltype(F)(0.5) * (KinF[I,J,2,iz+1] - KinF[I,J,1,iz+1])
    GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
    GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
    GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
    GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
    GradwF2 += -RhoCol * GradZ * dXdxI[3,3,2,ID,Iz,IF]
    GradwF1 += -RhoCol * GradZ * dXdxI[3,3,1,ID,Iz,IF]

    @atomic :monotonic F[Iz,ind,2] += (GraduF1 + GraduF2) / M[Iz,ind]
    @atomic :monotonic F[Iz,ind,3] += (GradvF1 + GradvF2) / M[Iz,ind]

    if Iz < Nz
      @atomic :monotonic F[Iz,ind,4] += GradwF2 / MRho[Iz,ind]
    end  
    if Iz > 1  
      @atomic :monotonic F[Iz-1,ind,4] += GradwF1 / MRho[Iz-1,ind]
    end  
  end

end  
