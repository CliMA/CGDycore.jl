@kernel inbounds = true function GradKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(X),@Const(M),@Const(Glob),Gravitation)


# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  pCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    pCol[I,J,iz] = p[Iz,ind]
    RhoCol[I,J,iz] = U[Iz,ind,1]
  end
  if iz == ColumnTilesDim && Iz < Nz
    pCol[I,J,iz+1] = p[Iz+1,ind]
    RhoCol[I,J,iz+1] = U[Iz+1,ind,1]
  end  

  @synchronize

  ID = I + (J - 1) * N  
  ind = Glob[ID,IF]

  if Iz <= Nz
    DXpCol = D[I,1] * pCol[1,J,iz]
    DYpCol = D[J,1] * pCol[I,1,iz]
    for k = 2 : N
      DXpCol += D[I,k] * pCol[k,J,iz]
      DYpCol += D[J,k] * pCol[I,k,iz]
    end
    @views Gradu, Gradv = Grad12(DXpCol,DYpCol,dXdxI[1:2,1:2,:,ID,Iz,IF]) 

    x = eltype(F)(0.5) * (X[ID,1,1,Iz,IF] + X[ID,2,1,Iz,IF])
    y = eltype(F)(0.5) * (X[ID,1,2,Iz,IF] + X[ID,2,2,Iz,IF])
    z = eltype(F)(0.5) * (X[ID,1,3,Iz,IF] + X[ID,2,3,Iz,IF])
    Grav = Gravitation(x,y,z)
    GradZ = -Grav * RhoCol[I,J,iz] *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,3,1,ID,Iz,IF] + dXdxI[3,3,2,ID,Iz,IF])
    Gradu += GradZ * (dXdxI[3,1,1,ID,Iz,IF] + dXdxI[3,1,2,ID,Iz,IF])
    Gradv += GradZ * (dXdxI[3,2,1,ID,Iz,IF] + dXdxI[3,2,2,ID,Iz,IF])
    @atomic :monotonic F[Iz,ind,2] += -Gradu / (M[Iz,ind,1] + M[Iz,ind,2]) / RhoCol[I,J,iz]
    @atomic :monotonic F[Iz,ind,3] += -Gradv / (M[Iz,ind,1] + M[Iz,ind,2]) / RhoCol[I,J,iz]
  end  

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
   
end

