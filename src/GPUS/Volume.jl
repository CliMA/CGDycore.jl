@kernel inbounds = true function VolumeNonLin!(F,@Const(U),@Const(p),@Const(D),@Const(DZ),
  @Const(dXdxI),@Const(JJ),@Const(X),@Const(Glob),CoriolisFun)

  K, I, J,  iz   = @index(Local, NTuple)
  _,_,_,Iz,IF = @index(Global, NTuple)

  M = @uniform @groupsize()[1]
  N = @uniform @groupsize()[2]
  Nz = @uniform @ndrange()[4]
  NF = @uniform @ndrange()[5]

  FCol = @localmem eltype(F) (M,N,N,3,5)
  DFCol = @localmem eltype(F) (M,N,N,5)

  iD = I + (J - 1) * N
  ind = Glob[iD,IF]
  if Iz <= Nz && IF <= NF
    @views FluxNonLin!(FCol[K,I,J,:,:],U[Iz,K,ind,:],p[Iz,K,ind],dXdxI[:,:,K,iD,Iz,IF])  
  end

  @synchronize
  iD = I + (J - 1) * N
  ind = Glob[iD,IF]
  x = X[iD,1,K,Iz,IF]
  y = X[iD,2,K,Iz,IF]
  z = X[iD,3,K,Iz,IF]
  FuCor, FvCor, FwCor = CoriolisFun(x,y,z,U[Iz,K,ind,2],U[Iz,K,ind,3],U[Iz,K,ind,4])
  for iv = 1 : 5
    f = D[I,1] * FCol[K,1,J,1,iv] + D[J,1] * FCol[K,I,1,2,iv]
    f1 = D[I,1] * FCol[K,1,J,1,iv] 
    f2 = D[J,1] * FCol[K,I,1,2,iv]
    for l = 2 : N
      f += D[I,l] * FCol[K,l,J,1,iv] + D[J,l] * FCol[K,I,l,2,iv]  
      f1 += D[I,l] * FCol[K,l,J,1,iv]
      f2 += D[J,l] * FCol[K,I,l,2,iv]  
    end
    for l = 1 : M
      f += DZ[K,l] * FCol[l,I,J,3,iv]  
    end
    if IF == 1 && iv == 2
      @show I,J,K,f,f1,f2,JJ[iD,K,Iz,IF]
    end  
    F[Iz,K,ind,iv] += f / JJ[iD,K,Iz,IF]
  end  
# F[Iz,K,ind,2] -= FuCor 
# F[Iz,K,ind,3] -= FvCor 
# F[Iz,K,ind,4] -= FwCor 
end

@inline function FluxNonLin!(F,U,p,dXdxI)
  RhouC = U[2]
  RhovC = U[3]
  RhowC = U[4]
  Rho = U[1]
  u = RhouC / Rho
  v = RhovC / Rho
  w = RhowC / Rho
  Th = U[5] / Rho

  F1 = -RhouC
  F2 = -RhovC
  F3 = -RhowC
  F[1,1] = dXdxI[1,1] * F1 + dXdxI[1,2] * F2 + dXdxI[1,3] * F3
  F[2,1] = dXdxI[2,1] * F1 + dXdxI[2,2] * F2 + dXdxI[2,3] * F3
  F[3,1] = dXdxI[3,1] * F1 + dXdxI[3,2] * F2 + dXdxI[3,3] * F3

  F1 = -RhouC * u - p
  F2 = -RhovC * u 
  F3 = -RhowC * u
  F[1,2] = dXdxI[1,1] * F1 + dXdxI[1,2] * F2 + dXdxI[1,3] * F3
  F[2,2] = dXdxI[2,1] * F1 + dXdxI[2,2] * F2 + dXdxI[2,3] * F3
  F[3,2] = dXdxI[3,1] * F1 + dXdxI[3,2] * F2 + dXdxI[3,3] * F3

  F1 = -RhouC * v 
  F2 = -RhovC * v -  p
  F3 = -RhowC * v
  F[1,3] = dXdxI[1,1] * F1 + dXdxI[1,2] * F2 + dXdxI[1,3] * F3
  F[2,3] = dXdxI[2,1] * F1 + dXdxI[2,2] * F2 + dXdxI[2,3] * F3
  F[3,3] = dXdxI[3,1] * F1 + dXdxI[3,2] * F2 + dXdxI[3,3] * F3

  F1 = -RhouC * w 
  F2 = -RhovC * w
  F3 = -RhowC * w - p
  F[1,4] = dXdxI[1,1] * F1 + dXdxI[1,2] * F2 + dXdxI[1,3] * F3
  F[2,4] = dXdxI[2,1] * F1 + dXdxI[2,2] * F2 + dXdxI[2,3] * F3
  F[3,4] = dXdxI[3,1] * F1 + dXdxI[3,2] * F2 + dXdxI[3,3] * F3


  F1 = -RhouC * Th
  F2 = -RhovC * Th
  F3 = -RhowC * Th
  F[1,5] = dXdxI[1,1] * F1 + dXdxI[1,2] * F2 + dXdxI[1,3] * F3
  F[2,5] = dXdxI[2,1] * F1 + dXdxI[2,2] * F2 + dXdxI[2,3] * F3
  F[3,5] = dXdxI[3,1] * F1 + dXdxI[3,2] * F2 + dXdxI[3,3] * F3



end


