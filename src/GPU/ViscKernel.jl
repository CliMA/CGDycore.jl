@kernel function HyperViscKernel!(F,MRho,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  if Iz <= Nz
    @views @inbounds uC, vC = Curl12(U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uCCol[I,J,iz] = uC
    @inbounds vCCol[I,J,iz] = vC
    @views @inbounds uD, vD = Contra12(U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uDCol[I,J,iz] = uD
    @inbounds vDCol[I,J,iz] = vD
    @inbounds ThCol[I,J,iz] = U[Iz,ind,5] / U[Iz,ind,1]
  end
  @synchronize

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @inbounds Dxc = D[I,1] * ThCol[1,J,iz]
    @inbounds Dyc = D[J,1] * ThCol[I,1,iz]
    @inbounds Curl[I,J,iz] = D[I,1] * uCCol[1,J,iz] + D[J,1] * vCCol[I,1,iz] 
    @inbounds Div[I,J,iz] = D[I,1] * uDCol[1,J,iz] + D[J,1] * vDCol[I,1,iz] 
    for k = 2 : N
      @inbounds Dxc += D[I,k] * ThCol[k,J,iz]
      @inbounds Dyc += D[J,k] * ThCol[I,k,iz] 
      @inbounds Curl[I,J,iz] += D[I,k] * uCCol[k,J,iz] + D[J,k] * vCCol[I,k,iz] 
      @inbounds Div[I,J,iz] += D[I,k] * uDCol[k,J,iz] + D[J,k] * vDCol[I,k,iz] 
    end
    @inbounds Curl[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @inbounds Div[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @views @inbounds (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views @inbounds (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds ThCxCol[I,J,iz] = tempx
    @inbounds ThCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @inbounds DxCurl = DW[I,1] * Curl[1,J,iz]
    @inbounds DyCurl = DW[J,1] * Curl[I,1,iz]
    @inbounds DxDiv = DW[I,1] * Div[1,J,iz]
    @inbounds DyDiv = DW[J,1] * Div[I,1,iz]
    @inbounds DivTh = DW[I,1] * ThCxCol[1,J,iz] + DW[J,1] * ThCyCol[I,1,iz]
    for k = 2 : N
      @inbounds DxCurl += DW[I,k] * Curl[k,J,iz]
      @inbounds DyCurl += DW[J,k] * Curl[I,k,iz]
      @inbounds DxDiv += DW[I,k] * Div[k,J,iz]
      @inbounds DyDiv += DW[J,k] * Div[I,k,iz]
      @inbounds DivTh += DW[I,k] * ThCxCol[k,J,iz] + DW[J,k] * ThCyCol[I,k,iz]
    end
    @views @inbounds FuC, FvC = Rot12(DxCurl,DyCurl,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @views @inbounds FuD, FvD = Grad12(DxDiv,DyDiv,dXdxI[1:2,1:2,:,ID,Iz,IF]) 
    @inbounds @atomic F[Iz,ind,1] += FuC / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,2] += FvC / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,3] += FuD / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,4] += FvD / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,5] += DivTh / M[Iz,ind]
    if Iz < Nz
      @inbounds @atomic MRho[Iz,ind] += U[Iz,ind,1] * JJ[ID,2,Iz,IF] 
    end  
    if Iz > 1
      @inbounds @atomic MRho[Iz-1,ind] += U[Iz,ind,1] * JJ[ID,1,Iz,IF] 
    end  
  end
end

@kernel function HyperViscKoeffKernel!(F,@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffCurl,KoeffGrad,KoeffDiv)

# gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @views @inbounds uC, vC = Curl12(Cache[Iz,ind,1],Cache[Iz,ind,2],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uCCol[I,J,iz] = uC
    @inbounds vCCol[I,J,iz] = vC
    @views @inbounds uD, vD = Contra12(Cache[Iz,ind,3],Cache[Iz,ind,4],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uDCol[I,J,iz] = uD
    @inbounds vDCol[I,J,iz] = vD
    @inbounds ThCol[I,J,iz] = Cache[Iz,ind,5]
  end
  @synchronize

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz
    @inbounds Dxc = D[I,1] * ThCol[1,J,iz]
    @inbounds Dyc = D[J,1] * ThCol[I,1,iz]
    @inbounds Curl[I,J,iz] = D[I,1] * uCCol[1,J,iz] + D[J,1] * vCCol[I,1,iz]
    @inbounds Div[I,J,iz] = D[I,1] * uDCol[1,J,iz] + D[J,1] * vDCol[I,1,iz]
    for k = 2 : N
      @inbounds Dxc += D[I,k] * ThCol[k,J,iz]
      @inbounds Dyc += D[J,k] * ThCol[I,k,iz]
      @inbounds Curl[I,J,iz] += D[I,k] * uCCol[k,J,iz] + D[J,k] * vCCol[I,k,iz]
      @inbounds Div[I,J,iz] += D[I,k] * uDCol[k,J,iz] + D[J,k] * vDCol[I,k,iz]
    end
    @inbounds Curl[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @inbounds Div[I,J,iz] /= (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @views @inbounds (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views @inbounds (tempx, tempy) = Contra12(U[Iz,ind,1],GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds ThCxCol[I,J,iz] = tempx
    @inbounds ThCyCol[I,J,iz] = tempy
  end

  @synchronize

  ID = I + (J - 1) * N
  @inbounds ind = Glob[ID,IF]
  if Iz <= Nz
    @inbounds DxCurl = DW[I,1] * Curl[1,J,iz]
    @inbounds DyCurl = DW[J,1] * Curl[I,1,iz]
    @inbounds DxDiv = DW[I,1] * Div[1,J,iz]
    @inbounds DyDiv = DW[J,1] * Div[I,1,iz]
    @inbounds DivTh = DW[I,1] * ThCxCol[1,J,iz] + DW[J,1] * ThCyCol[I,1,iz]
    for k = 2 : N
      @inbounds DxCurl += DW[I,k] * Curl[k,J,iz]
      @inbounds DyCurl += DW[J,k] * Curl[I,k,iz]
      @inbounds DxDiv += DW[I,k] * Div[k,J,iz]
      @inbounds DyDiv += DW[J,k] * Div[I,k,iz]
      @inbounds DivTh += DW[I,k] * ThCxCol[k,J,iz] + DW[J,k] * ThCyCol[I,k,iz]
    end
    @views @inbounds FuC, FvC = Rot12(DxCurl,DyCurl,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @views @inbounds FuD, FvD = Grad12(DxDiv,DyDiv,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds @atomic F[Iz,ind,2] += -(KoeffCurl * FuC + KoeffGrad * FuD) / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,3] += -(KoeffCurl * FvC + KoeffGrad * FvD) / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,5] += -KoeffDiv * DivTh / M[Iz,ind]
  end
end

@kernel function HyperViscTracerKernel!(FTr,@Const(Tr),@Const(Rho),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  TrCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCxCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCyCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  if Iz <= Nz && IF <= NF
    @inbounds TrCol[I,J,iz] = Tr[Iz,ind] / Rho[Iz,ind]
  end
  @synchronize

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz && IF <= NF
    @inbounds Dxc = D[I,1] * TrCol[1,J,iz]
    @inbounds Dyc = D[J,1] * TrCol[I,1,iz]
    for k = 2 : N
      @inbounds Dxc += D[I,k] * TrCol[k,J,iz]
      @inbounds Dyc += D[J,k] * TrCol[I,k,iz] 
    end
    @views @inbounds (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views @inbounds (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds TrCxCol[I,J,iz] = tempx
    @inbounds TrCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz <= Nz && IF <= NF
    @inbounds DivTr = DW[I,1] * TrCxCol[1,J,iz] + DW[J,1] * TrCyCol[I,1,iz]
    for k = 2 : N
      @inbounds DivTr += DW[I,k] * TrCxCol[k,J,iz] + DW[J,k] * TrCyCol[I,k,iz]
    end
    @inbounds @atomic FTr[Iz,ind] += DivTr / M[Iz,ind]
  end
end

@kernel function HyperViscWKernel!(Fw,@Const(w),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(MW),@Const(Glob)) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  wCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCxCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCyCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  if Iz < Nz && IF <= NF
    @inbounds wCol[I,J,iz] = w[Iz,ind] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz < Nz && IF <= NF
    @inbounds DxcF = D[I,1] * wCol[1,J,iz]
    @inbounds DycF = D[J,1] * wCol[I,1,iz]
    for k = 2 : N
      @inbounds DxcF += D[I,k] * wCol[k,J,iz]
      @inbounds DycF += D[J,k] * wCol[I,k,iz] 
    end
    @inbounds GradDx = ((dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    @inbounds GradDy = ((dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    @inbounds tempx = (dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * GradDy
    @inbounds tempy = (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * GradDy
    @inbounds wCxCol[I,J,iz] = tempx
    @inbounds wCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz < Nz && IF <= NF
    @inbounds Divw = DW[I,1] * wCxCol[1,J,iz] + DW[J,1] * wCyCol[I,1,iz]
    for k = 2 : N
      @inbounds Divw += DW[I,k] * wCxCol[k,J,iz] + DW[J,k] * wCyCol[I,k,iz]
    end
    @inbounds @atomic Fw[Iz,ind] += Divw / MW[Iz,ind]
  end
end

@kernel function HyperViscTracerKoeffKernel!(FTr,@Const(Cache),@Const(Rho),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffDiv) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  TrCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCxCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)
  TrCyCol = @localmem eltype(FTr) (N,N, ColumnTilesDim)

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz && IF <= NF
    @inbounds TrCol[I,J,iz] = Cache[Iz,ind] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz <= Nz && IF <= NF
    @inbounds Dxc = D[I,1] * TrCol[1,J,iz]
    @inbounds Dyc = D[J,1] * TrCol[I,1,iz]
    for k = 2 : N
      @inbounds Dxc += D[I,k] * TrCol[k,J,iz]
      @inbounds Dyc += D[J,k] * TrCol[I,k,iz] 
    end
    @views @inbounds (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views @inbounds (tempx, tempy) = Contra12(Rho[Iz,ind],GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds TrCxCol[I,J,iz] = tempx
    @inbounds TrCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz <= Nz && IF <= NF
    @inbounds DivTr = DW[I,1] * TrCxCol[1,J,iz] + DW[J,1] * TrCyCol[I,1,iz]
    for k = 2 : N
      @inbounds DivTr += DW[I,k] * TrCxCol[k,J,iz] + DW[J,k] * TrCyCol[I,k,iz]
    end
    @inbounds @atomic FTr[Iz,ind] += -KoeffDiv * DivTr / M[Iz,ind]
  end
end

@kernel function HyperViscWKoeffKernel!(Fw,@Const(w),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(MW),@Const(Glob),KoeffDivW) 

  I, J, iz   = @index(Local, NTuple)
  _,_,Iz,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  wCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCxCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  wCyCol = @localmem eltype(Fw) (N,N, ColumnTilesDim)
  if Iz < Nz && IF <= NF
    @inbounds wCol[I,J,iz] = w[Iz,ind] 
  end
  @synchronize

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]

  if Iz < Nz && IF <= NF
    @inbounds DxcF = D[I,1] * wCol[1,J,iz]
    @inbounds DycF = D[J,1] * wCol[I,1,iz]
    for k = 2 : N
      @inbounds DxcF += D[I,k] * wCol[k,J,iz]
      @inbounds DycF += D[J,k] * wCol[I,k,iz] 
    end
    @inbounds GradDx = ((dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    @inbounds GradDy = ((dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * DxcF +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * DycF) / (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    @inbounds tempx = (dXdxI[1,1,2,ID,Iz,IF] + dXdxI[1,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[1,2,2,ID,Iz,IF] + dXdxI[1,2,1,ID,Iz+1,IF]) * GradDy
    @inbounds tempy = (dXdxI[2,1,2,ID,Iz,IF] + dXdxI[2,1,1,ID,Iz+1,IF]) * GradDx +
      (dXdxI[2,2,2,ID,Iz,IF] + dXdxI[2,2,1,ID,Iz+1,IF]) * GradDy
    @inbounds wCxCol[I,J,iz] = tempx
    @inbounds wCyCol[I,J,iz] = tempy
  end

  @synchronize 

  ID = I + (J - 1) * N  
  @inbounds ind = Glob[ID,IF]
  if Iz < Nz && IF <= NF
    @inbounds Divw = DW[I,1] * wCxCol[1,J,iz] + DW[J,1] * wCyCol[I,1,iz]
    for k = 2 : N
      @inbounds Divw += DW[I,k] * wCxCol[k,J,iz] + DW[J,k] * wCyCol[I,k,iz]
    end
    @inbounds @atomic Fw[Iz,ind] += -KoeffDivW * Divw / MW[Iz,ind]
  end
end



















