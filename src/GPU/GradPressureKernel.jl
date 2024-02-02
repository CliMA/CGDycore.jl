@kernel function GradKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(MRho),@Const(Glob),Gravitation)

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  Pres = @localmem eltype(F) (N,N,2,ColumnTilesDim+2)

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @inbounds Pres[I,J,1,iz+1] = p[Iz,ind]
    @inbounds Pres[I,J,2,iz+1] = p[Iz,ind]
    if iz == 1 && Iz == 1
      Pres[I,J,2,iz] = p[Iz,ind] # oder Extrapolation ??   
    elseif iz == 1
      @inbounds Pres[I,J,2,iz] =  p[Iz-1,ind]
    end
    if iz == ColumnTilesDim && Iz < Nz
      @inbounds Pres[I,J,1,iz+2] = p[Iz+1,ind]
    end  
  end

  @synchronize


  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    x = X[ID,1,1,Iz,IF]
    y = X[ID,1,2,Iz,IF]
    z = X[ID,1,3,Iz,IF]
    Grav1 = Gravitation(x,y,z)
    x = X[ID,2,1,Iz,IF]
    y = X[ID,2,2,Iz,IF]
    z = X[ID,2,3,Iz,IF]
    Grav2 = Gravitation(x,y,z)
    @inbounds DXPres1 = D[I,1] * Pres[1,J,1,iz+1]
    @inbounds DYPres1 = D[J,1] * Pres[I,1,1,iz+1]
    @inbounds DXPres2 = D[I,1] * Pres[1,J,2,iz+1]
    @inbounds DYPres2 = D[J,1] * Pres[I,1,2,iz+1]
    for k = 2 : N
      @inbounds DXPres1 += D[I,k] * Pres[k,J,1,iz+1]
      @inbounds DYPres1 += D[J,k] * Pres[I,k,1,iz+1]
      @inbounds DXPres2 += D[I,k] * Pres[k,J,2,iz+1]
      @inbounds DYPres2 += D[J,k] * Pres[I,k,2,iz+1]
    end
    @inbounds GraduF1 =
      -(dXdxI[1,1,1,ID,Iz,IF]  * DXPres1 + dXdxI[2,1,1,ID,Iz,IF]  * DYPres1)
    @inbounds GradvF1 =
      -(dXdxI[1,2,1,ID,Iz,IF]  * DXPres1 + dXdxI[2,2,1,ID,Iz,IF]  * DYPres1)
    @inbounds GraduF2 =
      -(dXdxI[1,1,2,ID,Iz,IF]  * DXPres2 + dXdxI[2,1,2,ID,Iz,IF]  * DYPres2)
    @inbounds GradvF2 =
      -(dXdxI[1,2,2,ID,Iz,IF]  * DXPres2 + dXdxI[2,2,2,ID,Iz,IF]  * DYPres2)

    GradwF1 = eltype(F)(0)
    GradwF2 = eltype(F)(0)
    RhoCol = U[Iz,ind,1]
    if Iz > 1 
      @inbounds GradZ = -Grav1 * RhoCol * JJ[ID,1,Iz,IF] / dXdxI[3,3,1,ID,Iz,IF]
      @inbounds GraduF1 += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
      @inbounds GradvF1 += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
      @inbounds GradZ = eltype(F)(0.5) * (Pres[I,J,1,iz+1] - Pres[I,J,2,iz] )  
      @inbounds GradwF1 += -GradZ * dXdxI[3,3,1,ID,Iz,IF]
    end  
    if Iz < Nz
      @inbounds GradZ = -Grav2 * RhoCol * JJ[ID,2,Iz,IF] / dXdxI[3,3,2,ID,Iz,IF]
      @inbounds GraduF2 += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
      @inbounds GradvF2 += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
      @inbounds GradZ = eltype(F)(0.5) * (Pres[I,J,1,iz+2] - Pres[I,J,2,iz+1] )  
      @inbounds GradwF2 += -GradZ * dXdxI[3,3,2,ID,Iz,IF]
    end  

    @inbounds @atomic F[Iz,ind,2] += (GraduF1 + GraduF2) / M[Iz,ind] / RhoCol
    @inbounds @atomic F[Iz,ind,3] += (GradvF1 + GradvF2) / M[Iz,ind] / RhoCol

    if Iz < Nz
      @inbounds @atomic F[Iz,ind,4] += (GradwF2 - Grav2 * RhoCol * JJ[ID,2,Iz,IF])/ MRho[Iz,ind]
    end  
    if Iz > 1  
      @inbounds @atomic F[Iz-1,ind,4] += (GradwF1 - Grav1 * RhoCol * JJ[ID,1,Iz,IF]) / MRho[Iz-1,ind]
    end  
  end

end  
