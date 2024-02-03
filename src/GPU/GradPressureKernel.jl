@kernel function GradKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(MRho),@Const(Glob),Gravitation)


# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  Pres = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @inbounds Pres[I,J,iz] = p[Iz,ind]
  end
  if iz == ColumnTilesDim && Iz < Nz
    @inbounds Pres[I,J,iz+1] = p[Iz+1,ind]
  end  

  @synchronize

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @inbounds DXPres = D[I,1] * Pres[1,J,iz]
    @inbounds DYPres = D[J,1] * Pres[I,1,iz]
    for k = 2 : N
      @inbounds DXPres += D[I,k] * Pres[k,J,iz]
      @inbounds DYPres += D[J,k] * Pres[I,k,iz]
    end
    @views @inbounds Gradu, Gradv = Grad12(DXPres,DYPres,dXdxI[1:2,1:2,:,ID,Iz,IF]) 

    x = eltype(F)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(F)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(F)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    Grav = Gravitation(x,y,z)
    @inbounds GradZ = -Grav * U[Iz,ind,1] *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,3,1,ID,Iz,IF] + dXdxI[3,3,2,ID,Iz,IF])
    @inbounds Gradu += GradZ * (dXdxI[3,1,1,ID,Iz,IF] + dXdxI[3,1,2,ID,Iz,IF])
    @inbounds Gradv += GradZ * (dXdxI[3,2,1,ID,Iz,IF] + dXdxI[3,2,2,ID,Iz,IF])
    @inbounds @atomic F[Iz,ind,2] += -Gradu / M[Iz,ind] / U[Iz,ind,1]
    @inbounds @atomic F[Iz,ind,3] += -Gradv / M[Iz,ind] / U[Iz,ind,1]
  end  

  if Iz < Nz
    @inbounds GradZ = eltype(F)(0.5) * (Pres[I,J,iz+1] - Pres[I,J,iz])  
    @inbounds Gradw =  GradZ* (dXdxI[3,3,2,ID,Iz,IF] + dXdxI[3,3,1,ID,Iz+1,IF])
    x = eltype(F)(0.5) * (X[ID,2,1,Iz,IF] + X[ID,2,1,Iz+1,IF])
    y = eltype(F)(0.5) * (X[ID,2,2,Iz,IF] + X[ID,1,2,Iz+1,IF])
    z = eltype(F)(0.5) * (X[ID,2,3,Iz,IF] + X[ID,1,3,Iz+1,IF])
    Grav = Gravitation(x,y,z)
    @inbounds @atomic F[Iz,ind,4] += -(Gradw + 
      Grav * (U[Iz,ind,1] * JJ[ID,2,Iz,IF] + U[Iz+1,ind,1] * JJ[ID,1,Iz+1,IF])) /
      MRho[Iz,ind]
  end      
   
end

