@kernel function MomentumKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(MRho),@Const(M),@Const(Glob),Phys)

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  
  RhoCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  uCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds RhoCol[I,J,iz] = U[Iz,ind,1]
    @inbounds uCol[I,J,iz] = U[Iz,ind,2]
    @inbounds vCol[I,J,iz] = U[Iz,ind,3]
    @inbounds wCol[I,J,iz+1] = U[Iz,ind,4]
    if Iz == 1
      wCol[I,J,1] = -(dXdxI[3,1,1,ID,1,IF] * U[Iz,ind,2] + 
        dXdxI[3,2,1,ID,1,IF] * U[Iz,ind,3]) / dXdxI[3,3,1,ID,1,IF]
     elseif iz == 1
       wCol[I,J,1] = U[Iz-1,ind,4] 
    end    
  end  

  @synchronize
    
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds uCon1 = -RhoCol[I,J,iz] * (dXdxI[1,1,1,ID,Iz,IF] * uCol[I,J,iz] +
      dXdxI[1,2,1,ID,Iz,IF] * vCol[I,J,iz] + dXdxI[1,3,1,ID,Iz,IF] * wCol[I,J,iz])
    @inbounds uCon2 = -RhoCol[I,J,iz] * (dXdxI[1,1,2,ID,Iz,IF] * uCol[I,J,iz] +
      dXdxI[1,2,2,ID,Iz,IF] * vCol[I,J,iz] + dXdxI[1,3,2,ID,Iz,IF] * wCol[I,J,iz+1])
    @inbounds vCon1 = -RhoCol[I,J,iz] * (dXdxI[2,1,1,ID,Iz,IF] * uCol[I,J,iz] +
      dXdxI[2,2,1,ID,Iz,IF] * vCol[I,J,iz] + dXdxI[2,3,1,ID,Iz,IF] * wCol[I,J,iz])
    @inbounds vCon2 = -RhoCol[I,J,iz] * (dXdxI[2,1,2,ID,Iz,IF] * uCol[I,J,iz] +
      dXdxI[2,2,2,ID,Iz,IF] * vCol[I,J,iz] + dXdxI[2,3,2,ID,Iz,IF] * wCol[I,J,iz+1])
    @inbounds wCon1 = -RhoCol[I,J,iz] * (dXdxI[3,1,1,ID,Iz,IF] * uCol[I,J,iz] +
      dXdxI[3,2,1,ID,Iz,IF] * vCol[I,J,iz] + dXdxI[3,3,1,ID,Iz,IF] * wCol[I,J,iz])
    @inbounds wCon2 = -RhoCol[I,J,iz] * (dXdxI[3,1,2,ID,Iz,IF] * uCol[I,J,iz] +
      dXdxI[3,2,2,ID,Iz,IF] * vCol[I,J,iz] + dXdxI[3,3,2,ID,Iz,IF] * wCol[I,J,iz+1])

    @inbounds Dxu = D[I,1] * uCol[1,J,iz]
    @inbounds Dyu = D[J,1] * uCol[I,1,iz]
    @inbounds Dxv = D[I,1] * vCol[1,J,iz]
    @inbounds Dyv = D[J,1] * vCol[I,1,iz]
    @inbounds Dxw1 = D[I,1] * wCol[1,J,iz]
    @inbounds Dyw1 = D[J,1] * wCol[I,1,iz]
    @inbounds Dxw2 = D[I,1] * wCol[1,J,iz+1]
    @inbounds Dyw2 = D[J,1] * wCol[I,1,iz+1]
    Izp = min(Iz+1,Nz)
    Izm = max(Iz-1,1)
    ind = Glob[ID,IF]
    Dzu2 = 1/2 * (U[Izp,ind,2] - uCol[I,J,iz])
    Dzv2 = 1/2 * (U[Izp,ind,3] - vCol[I,J,iz])
    Dzu1 = 1/2 * (uCol[I,J,iz] - U[Izm,ind,2])
    Dzv1 = 1/2 * (vCol[I,J,iz] - U[Izm,ind,3])
    Dzw = 1/2 * (wCol[I,J,iz+1] - wCol[I,J,iz]) 
    for k = 2 : N
      @inbounds Dxu += D[I,k] * uCol[k,J,iz]
      @inbounds Dyu += D[J,k] * uCol[I,k,iz]
      @inbounds Dxv += D[I,k] * vCol[k,J,iz]
      @inbounds Dyv += D[J,k] * vCol[I,k,iz]
      @inbounds Dxw1 += D[I,k] * wCol[k,J,iz]
      @inbounds Dyw1 += D[J,k] * wCol[I,k,iz]
      @inbounds Dxw2 += D[I,k] * wCol[k,J,iz+1]
      @inbounds Dyw2 += D[J,k] * wCol[I,k,iz+1]
    end  

    @inbounds @atomic F[Iz,ind,2] += ((uCon1 + uCon2) * Dxu + (vCon1 + vCon2) * Dyu + 
    wCon1 * Dzu1 + wCon2 * Dzu2) / M[Iz,ind] / RhoCol[I,J,iz]
    @inbounds @atomic F[Iz,ind,3] += ((uCon1 + uCon2) * Dxv + (vCon1 + vCon2) * Dyv +
    wCon1 * Dzv1 + wCon2 * Dzv2) / M[Iz,ind] / RhoCol[I,J,iz]
  end  
  if Iz > 1
    @inbounds @atomic F[Iz-1,ind,4] += (uCon1 * Dxw1 + vCon1 * Dyw1 + wCon1 * Dzw) / MRho[Iz-1,ind] 
  end  
  if Iz < Nz
    @inbounds @atomic F[Iz,ind,4] += (uCon2 * Dxw2 + vCon2 * Dyw2 + wCon2 * Dzw) / MRho[Iz,ind]
  end  
end  

@kernel function GradKernel!(F,@Const(U),@Const(p),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(MRho),@Const(Glob),Phys,::Val{BANK}=Val(1)) where BANK

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  Pres = @localmem eltype(F) (N,N,ColumnTilesDim+1)

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds Pres[I,J,iz] = p[Iz,ind]
  end
  if iz == ColumnTilesDim && Iz < Nz
    @inbounds Pres[I,J,iz+1] = p[Iz+1,ind]
  end  

  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds DXPres = D[I,1] * Pres[1,J,iz]
    @inbounds DYPres = D[J,1] * Pres[I,1,iz]
    for k = 2 : N
      @inbounds DXPres += D[I,k] * Pres[k,J,iz]
      @inbounds DYPres += D[J,k] * Pres[I,k,iz]
    end
    @views @inbounds Gradu, Gradv = Grad12(DXPres,DYPres,dXdxI[1:2,1:2,:,ID,Iz,IF]) 
    @views @inbounds Gradw1, Gradw2 = Grad3(DXPres,DYPres,dXdxI[1:3,1:3,:,ID,Iz,IF]) 

    @inbounds ind = Glob[ID,IF]
    @inbounds GradZ = -Phys.Grav * U[Iz,ind,1] *
        (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]) / (dXdxI[3,3,1,ID,Iz,IF] + dXdxI[3,3,2,ID,Iz,IF])
    @inbounds Gradu += GradZ * (dXdxI[3,1,1,ID,Iz,IF] + dXdxI[3,1,2,ID,Iz,IF])
    @inbounds Gradv += GradZ * (dXdxI[3,2,1,ID,Iz,IF] + dXdxI[3,2,2,ID,Iz,IF])
    @inbounds @atomic F[Iz,ind,2] += -Gradu / M[Iz,ind] / U[Iz,ind,1]
    @inbounds @atomic F[Iz,ind,3] += -Gradv / M[Iz,ind] / U[Iz,ind,1]
    if Iz > 1
      @inbounds @atomic F[Iz-1,ind,4] += -Gradw1 / MRho[Iz-1,ind]
    end  
  end  

  if Iz < Nz
    @inbounds GradZ = 1/2 * (Pres[I,J,iz+1] - Pres[I,J,iz])  
    @inbounds Gradw =  GradZ* (dXdxI[3,3,2,ID,Iz,IF] + dXdxI[3,3,1,ID,Iz+1])
    @inbounds @atomic F[Iz,ind,4] += -(Gradw + Gradw2 +
      Phys.Grav * (U[Iz,ind,1] * JJ[ID,2,Iz,IF] + U[Iz+1,ind,1] * JJ[ID,1,Iz+1,IF])) /
      MRho[Iz,ind]
  end      
   
end

@kernel function RhoGradKinKernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N,ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N,ColumnTilesDim+1)
  GraduF = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  GradvF = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  GradwF = @localmem eltype(F) (N,N,2,ColumnTilesDim)
  KinF = @localmem eltype(F) (N,N,2,ColumnTilesDim)

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds RhoCol[I,J,iz] = U[Iz,ind,1]
    @inbounds uCol[I,J,iz] = U[Iz,ind,2]
    @inbounds vCol[I,J,iz] = U[Iz,ind,3]
    @inbounds wCol[I,J,iz+1] = U[Iz,ind,4]
    @inbounds @views @. GraduF[I,J,:,iz] = 0
    @inbounds @views @. GradvF[I,J,:,iz] = 0
    @inbounds @views @. GradwF[I,J,:,iz] = 0
  end

  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    KinF[I,J,1,iz] = 1/2 * (uCol[I,J,iz] * uCol[I,J,iz] + vCol[I,J,iz] * vCol[I,J,iz])  
    KinF[I,J,2,iz] = KinF[I,J,1,iz] + 1/2 * wCol[I,J,iz+1] * wCol[I,J,iz+1]
    KinF[I,J,1,iz] +=  1/2 * wCol[I,J,iz] * wCol[I,J,iz]
  end  

  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz

  DXKinF1 = 0
  DYKinF1 = 0
  DXKinF2 = 0
  DYKinF2 = 0
  I = mod(ID-1,N) + 1
  J = div(ID-I,N) + 1
  for k = 1 : N
    @inbounds DXKinF1 += D[I,k] * KinF[k,J,1,iz]
    @inbounds DYKinF1 += D[J,k] * KinF[I,k,1,iz]
    @inbounds DXKinF2 += D[I,k] * KinF[k,J,2,iz]
    @inbounds DYKinF2 += D[J,k] * KinF[I,k,2,iz]
  end
  @inbounds GraduF[I,J,1,iz] +=
      -RhoCol[I,J,iz] * (dXdxI[1,1,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,1,1,ID,Iz,IF]  * DYKinF1)
  @inbounds GradvF[I,J,1,iz] +=
      -RhoCol[I,J,iz] * (dXdxI[1,2,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,2,1,ID,Iz,IF]  * DYKinF1)
  @inbounds GradwF[I,J,1,iz] +=
      -RhoCol[I,J,iz] * (dXdxI[1,3,1,ID,Iz,IF]  * DXKinF1 + dXdxI[2,3,1,ID,Iz,IF]  * DYKinF1)
  @inbounds GraduF[I,J,2,iz] +=
      -RhoCol[I,J,iz] * (dXdxI[1,1,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,1,2,ID,Iz,IF]  * DYKinF2)
  @inbounds GradvF[I,J,2,iz] +=
      -RhoCol[I,J,iz] * (dXdxI[1,2,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,2,2,ID,Iz,IF]  * DYKinF2)
  @inbounds GradwF[I,J,2,iz] +=
      -RhoCol[I,J,iz] * (dXdxI[1,3,2,ID,Iz,IF]  * DXKinF2 + dXdxI[2,3,1,ID,Iz,IF]  * DYKinF2)
  if iz > 1    
    @inbounds GraduZ11 = 1/2 * KinF[I,J,1,iz] * dXdxI[3,1,2,ID,Iz-1,IF]
    @inbounds GradvZ11 = 1/2 * KinF[I,J,1,iz] * dXdxI[3,2,2,ID,Iz-1,IF]
    @inbounds @atomic GraduF[I,J,2,iz-1] += -RhoCol[I,J,iz-1] * GraduZ11
    @inbounds @atomic GradvF[I,J,2,iz-1] += -RhoCol[I,J,iz-1] * GradvZ11
  end  
  @inbounds GraduZ12 = 1/2 * KinF[I,J,1,iz] * dXdxI[3,1,1,ID,Iz,IF]
  @inbounds GradvZ12 = 1/2 * KinF[I,J,1,iz] * dXdxI[3,2,1,ID,Iz,IF]
  @inbounds GraduZ21 = 1/2 * KinF[I,J,2,iz] * dXdxI[3,1,2,ID,Iz,IF]
  @inbounds GradvZ21 = 1/2 * KinF[I,J,2,iz] * dXdxI[3,2,2,ID,Iz,IF]
  @inbounds @atomic GraduF[I,J,1,iz] += RhoCol[I,J,iz] * GraduZ12 
  @inbounds @atomic GraduF[I,J,1,iz] += -RhoCol[I,J,iz] * GraduZ21 
  @inbounds @atomic GradvF[I,J,1,iz] += RhoCol[I,J,iz] * GradvZ12 
  @inbounds @atomic GradvF[I,J,1,iz] += -RhoCol[I,J,iz] * GradvZ21 
  if Iz < Nz
    @inbounds GraduZ22 = 1/2 * KinF[I,J,2,iz] * dXdxI[3,1,1,ID,Iz+1,IF]
    @inbounds GradvZ22 = 1/2 * KinF[I,J,2,iz] * dXdxI[3,2,1,ID,Iz+1,IF]
    @inbounds @atomic GraduF[I,J,2,iz+1] += RhoCol[I,J,iz+1] * GraduZ22
    @inbounds @atomic GradvF[I,J,2,iz+1] += RhoCol[I,J,iz+1] * GradvZ22
  end  

  @inbounds GradZ = 1/2 * RhoCol[I,J,iz] * (KinF[I,J,2,iz] - KinF[I,J,1,iz])
  @inbounds GraduF[I,J,2,iz] += -GradZ * dXdxI[3,1,2,ID,Iz,IF]
  @inbounds GraduF[I,J,1,iz] += -GradZ * dXdxI[3,1,1,ID,Iz,IF]
  @inbounds GradvF[I,J,2,iz] += -GradZ * dXdxI[3,2,2,ID,Iz,IF]
  @inbounds GradvF[I,J,1,iz] += -GradZ * dXdxI[3,2,1,ID,Iz,IF]
  @inbounds GradwF[I,J,2,iz] += -GradZ * dXdxI[3,3,2,ID,Iz,IF]
  @inbounds GradwF[I,J,1,iz] += -GradZ * dXdxI[3,3,1,ID,Iz,IF]

  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind,2] += (GraduF[I,J,1,iz] + GraduF[I,J,2,iz]) / M[Iz,ind] / U[Iz,ind,1]
    @inbounds @atomic F[Iz,ind,3] += (GraduF[I,J,1,iz] + GraduF[I,J,2,iz]) / M[Iz,ind] / U[Iz,ind,1]
    if iz > 1
      @inbounds @atomic F[Iz,ind,4] += (GradwF[I,J,2,iz-1] + GradwF[I,J,1,iz]) / 
        (M[Iz,ind] * U[Iz,ind,1] + M[Iz-1,ind] * U[Iz-1,ind,1])
    end  
  end
end

@kernel function DivRhoGradKernel!(F,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),::Val{BANK}=Val(1)) where BANK

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim)
  FCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    @inbounds ind = Glob[ID,IF]
    @inbounds cCol[I,J,iz] = U[Iz,ind,5] / U[Iz,ind,1]
    @inbounds FCol[I,J,iz] = 0.0
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    Dxc = D[I,1] * cCol[1,J,iz]
    Dyc = D[J,1] * cCol[I,1,iz]
    for k = 2 : N
      @inbounds Dxc = Dxc + D[I,k] * cCol[k,J,iz]
      @inbounds Dyc = Dyc + D[J,k] * cCol[I,k,iz] 
    end
    @views @inbounds (GradDx, GradDy) = Grad12(Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views @inbounds (tempx, tempy) = Contra12(GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @inbounds @atomic FCol[k,J,iz] += DW[k,I] * tempx
      @inbounds @atomic FCol[I,k,iz] += DW[k,J] * tempy
    end
  end

  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    @inbounds ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind,5] += FCol[I,J,iz] / M[Iz,ind]
  end
end

@kernel function HyperViscKernel!(F,MRho,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob)) 

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @views @inbounds uC, vC = Curl12(U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uCCol[I,J,iz] = uC
    @inbounds vCCol[I,J,iz] = vC
    @views @inbounds uD, vD = Contra12(U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uDCol[I,J,iz] = uD
    @inbounds vDCol[I,J,iz] = vD
    @inbounds ThCol[I,J,iz] = U[Iz,ind,5] / U[Iz,ind,1]
  end
  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
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
  Iz = (gz - 1) * ColumnTilesDim + iz
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
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
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

@kernel function HyperViscKernelKoeff!(F,@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),KoeffCurl,KoeffGrad,KoeffDiv) 

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vDCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Curl = @localmem eltype(F) (N,N, ColumnTilesDim)
  Div = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCxCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCyCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @views @inbounds uC, vC = Curl12(Cache[Iz,ind,1],Cache[Iz,ind,2],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uCCol[I,J,iz] = uC
    @inbounds vCCol[I,J,iz] = vC
    @views @inbounds uD, vD = Contra12(Cache[Iz,ind,3],Cache[Iz,ind,4],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uDCol[I,J,iz] = uD
    @inbounds vDCol[I,J,iz] = vD
    @inbounds ThCol[I,J,iz] = Cache[Iz,ind,5] 
  end
  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
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
    @inbounds ind = Glob[ID,IF]
    @views @inbounds (tempx, tempy) = Contra12(U[Iz,ind,1],GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds ThCxCol[I,J,iz] = tempx
    @inbounds ThCyCol[I,J,iz] = tempy
  end

  @synchronize 
  Iz = (gz - 1) * ColumnTilesDim + iz
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
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @views @inbounds FuC, FvC = Rot12(DxCurl,DyCurl,dXdxI[1:2,1:2,:,ID,Iz,IF])
    @views @inbounds FuD, FvD = Grad12(DxDiv,DyDiv,dXdxI[1:2,1:2,:,ID,Iz,IF]) 
    @inbounds @atomic F[Iz,ind,2] += -(KoeffCurl * FuC + KoeffGrad * FuD) / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,3] += -(KoeffCurl * FvC + KoeffGrad * FvD) / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,5] += -KoeffDiv * DivTh / M[Iz,ind]
  end
end

@kernel function DivRhoGradKernel1!(F,@Const(U),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  ThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FRhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FuCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FvCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FThCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  CurlCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    @inbounds ind = Glob[ID,IF]
    @inbounds RhoCol[I,J,iz+1] = U[Iz,ind,1]
    @inbounds uCol[I,J,iz+1] = U[Iz,ind,2]
    @inbounds vCol[I,J,iz+1] = U[Iz,ind,3]
    @inbounds wCol[I,J,iz+1] = U[Iz,ind,4]
    @inbounds ThCol[I,J,iz+1] = U[Iz,ind,5] / RhoCol[I,J,iz+1]
    @inbounds FRhoCol[I,J,iz+1] = 0
    @inbounds FThCol[I,J,iz+1] = 0
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
#   DivGrad Th
    ID = I + (J - 1) * N  
    Dxc = 0
    Dyc = 0
    for k = 1 : N
      @inbounds Dxc = Dxc + D[I,k] * ThCol[k,J,iz]
      @inbounds Dyc = Dyc + D[J,k] * ThCol[I,k,iz] 
    end
    @inbounds GradDx = ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * Dxc +
      (dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * Dyc) / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @inbounds GradDy = ((dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * Dxc +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * Dyc) / (JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF])
    @inbounds tempx = (dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * GradDx +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * GradDy
    @inbounds tempy = (dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * GradDx +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * GradDy
    for k = 1 : N
      @inbounds @atomic FThCol[k,J,iz] += DW[k,I] * tempx
      @inbounds @atomic FThCol[I,k,iz] += DW[k,J] * tempy
    end
#   Curl (u,v)
    @inbounds tempx = (dXdxI[I,J,1,iz,1,1,IF] + dXdxI[I,J,2,iz,1,1],IF) * vC[I,J,iz] -
      (dXdxI[I,J,1,iz,1,2,IF] + dXdxI[I,J,2,iz,1,2,IF]) * uC[I,J,iz]
    @views @. tempy = (dXdxI[I,J,1,iz,2,1,IF] + dXdxI[I,J,2,iz,2,1,IF]) * vC[I,J,iz] -
      (dXdxI[I,J,1,iz,2,2,IF] + dXdxI[I,J,2,iz,2,2,IF]) * uC[I,J,iz]
    for k = 1 : N
      @inbounds @atomic CurlCol[k,J,iz] += D[k,I] * tempx + D[k,J] * tempy
    end
  end

  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
#   CurlCurl (u,v)
    I = mod(ID-1,N) + 1
    J = div(ID-I,N) + 1
    CurlCol[I,J,iz] /= (J[I,J,1,iz] + J[I,J,2,iz])
    DxCurl = eltype(F)(0)
    DyCurl = eltype(F)(0)
    for k = 1 : N
      @inbounds DxCurl += DW[I,k] * CurlCol[k,J,iz]
      @inbounds DyCurl += DW[J,k] * CurlCol[I,k,iz] 
    end

    @inbounds FvCol[I,J,iz] = (-(dXdxI[I,J,1,iz,1,1,IF] + dXdxI[I,J,2,iz,1,1,IF]) * DxCurl -
      (dXdxI[I,J,1,iz,2,1,IF] + dXdxI[I,J,2,iz,2,1,IF]) * DyCurl)
    @inbounds FuCol[I,J,iz] = ((dXdxI[I,J,1,iz,1,2,IF] + dXdxI[I,J,2,iz,1,2,IF]) * DxCurl +
      (dXdxI[I,J,1,iz,2,2,IF] + dXdxI[I,J,2,iz,2,2,IF]) * DyCurl)
  end

  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind,5] += FThCol[I,J,iz] / M[Iz,ind]
  end
end

@kernel function DivRhoTrCentralKernel!(F,@Const(c),@Const(uC),@Const(vC),@Const(w),
  @Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),::Val{BANK}=Val(1)) where BANK

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  uCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  vCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  wCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+1)
  FCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds wCol[I,J,iz+1] = w[Iz,ind]
    @inbounds cCol[I,J,iz+1] = c[Iz,ind]
    @inbounds uCol[I,J,iz+1] = uC[Iz,ind]
    @inbounds vCol[I,J,iz+1] = vC[Iz,ind]
    @inbounds FCol[I,J,iz+1] = 0
    if iz == 1 && Iz > 1
      @inbounds cCol[I,J,1] = c[Iz-1,ind]
      @inbounds uCol[I,J,1] = uC[Iz-1,ind]
      @inbounds vCol[I,J,1] = vC[Iz-1,ind]
      @inbounds wCol[I,J,1] = w[Iz,ind]
      @inbounds FCol[I,J,1] = 0
    elseif iz == 1 && Iz == 1
      @inbounds cCol[I,J,1] = c[1,ind]
      @inbounds wCol[I,J,1] = 0
      @inbounds FCol[I,J,1] = 0
    end
    if iz == ColumnTilesDim && Iz < Nz
      @inbounds cCol[I,J,ColumnTilesDim+2] = c[Iz+1,ind]
      @inbounds uCol[I,J,ColumnTilesDim+2] = uC[Iz+1,ind]
      @inbounds vCol[I,J,ColumnTilesDim+2] = vC[Iz+1,ind]
      @inbounds FCol[I,J,ColumnTilesDim+2] = 0
    elseif iz == ColumnTilesDim && Iz == Nz
      @inbounds cCol[I,J,ColumnTilesDim+2] = c[Nz,ind]
      @inbounds FCol[I,J,ColumnTilesDim+2] = 0
    end
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz < Nz 
    @inbounds wCon = dXdxI[I,J,2,Iz,3,1,IF] * uCol[I,J,iz+1] + 
      dXdxI[I,J,2,Iz,3,2,IF] * vCol[I,J,iz+1] + 
      dXdxI[I,J,1,Iz+1,3,1,IF] * uCol[I,J,iz+2] + 
      dXdxI[I,J,2,Iz+1,3,2,IF] * vCol[I,J,iz+2] + 
      (dXdxI[I,J,2,Iz,3,3,IF] + dXdxI[I,J,1,Iz+1,3,3,IF]) * wCol[I,J,iz+1]
    @inbounds cF = (JJ[ID,2,Iz,IF] * cCol[I,J,iz+1] + JJ[ID,1,Iz+1,IF] * cCol[I,J,iz+2]) /
      (JJ[ID,2,Iz,IF] + JJ[ID,1,Iz+1,IF])
    Flux = 0.5 * wCon * cF
    @inbounds @atomic FCol[I,J,iz+1] += -Flux
    @inbounds @atomic FCol[I,J,iz+2] += Flux
  end 

  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds tempx = -cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,1,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,1,3,IF] * wCol[I,J,iz+1])
    @inbounds tempy = -cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,2,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,2,3,IF] * wCol[I,J,iz+1])
    for k = 1 : N
      @inbounds @atomic FCol[k,J,iz+1] += D[k,I] * tempx
      @inbounds @atomic FCol[I,k,iz+1] += D[k,J] * tempy
    end
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz 
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind] += FCol[I,J,iz+1] / M[Iz,ind]
    if iz == 1 && Iz >  1
      @inbounds @atomic F[Iz-1,ind] += FCol[I,J,iz] / M[Iz-1,ind]
    end
    if iz == ColumnTilesDim && Iz <  Nz
      @inbounds @atomic F[Iz+1,ind] += FCol[I,J,iz+2] / M[Iz+1,ind]
    end
  end
end

@kernel function DivRhoTrUpwindKernel!(F,@Const(c),@Const(Rho),@Const(uC),@Const(vC),@Const(w),
  @Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),::Val{BANK}=Val(1)) where BANK

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  cCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  uCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  vCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  RhoCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  wCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+1)
  FCol = @localmem eltype(F) (N+BANK,N, ColumnTilesDim+2)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds wCol[I,J,iz+1] = w[Iz,ind]
    @inbounds RhoCol[I,J,iz+1] = Rho[Iz,ind]
    @inbounds cCol[I,J,iz+1] = c[Iz,ind] / RhoCol[I,J,iz+1]
    @inbounds uCol[I,J,iz+1] = uC[Iz,ind]
    @inbounds vCol[I,J,iz+1] = vC[Iz,ind]
    @inbounds FCol[I,J,iz+1] = 0
    if iz == 1 && Iz > 1
      @inbounds RhoCol[I,J,1] = Rho[Iz-1,ind]
      @inbounds cCol[I,J,1] = c[Iz-1,ind] / RhoCol[I,J,1]
      @inbounds uCol[I,J,1] = uC[Iz-1,ind]
      @inbounds vCol[I,J,1] = vC[Iz-1,ind]
      @inbounds wCol[I,J,1] = w[Iz,ind]
      @inbounds FCol[I,J,1] = 0
    elseif iz == 1 && Iz == 1
      @inbounds RhoCol[I,J,1] = Rho[1,ind]
      @inbounds cCol[I,J,1] = c[1,ind] / RhoCol[I,J,1]
      @inbounds wCol[I,J,1] = 0
      @inbounds FCol[I,J,1] = 0
    end
    if iz == ColumnTilesDim && Iz < Nz
      @inbounds RhoCol[I,J,ColumnTilesDim+2] = Rho[Iz+1,ind]
      @inbounds cCol[I,J,ColumnTilesDim+2] = c[Iz+1,ind] / RhoCol[I,J,ColumnTilesDim+2]
      @inbounds uCol[I,J,ColumnTilesDim+2] = uC[Iz+1,ind]
      @inbounds vCol[I,J,ColumnTilesDim+2] = vC[Iz+1,ind]
      @inbounds FCol[I,J,ColumnTilesDim+2] = 0
    elseif iz == ColumnTilesDim && Iz == Nz
      @inbounds RhoCol[I,J,ColumnTilesDim+2] = Rho[Nz,ind]
      @inbounds cCol[I,J,ColumnTilesDim+2] = c[Nz,ind] / RhoCol[I,J,ColumnTilesDim+2]
      @inbounds FCol[I,J,ColumnTilesDim+2] = 0
    end
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz < Nz 
    @inbounds wCon = RhoCol[I,J,iz+1] * (dXdxI[I,J,2,Iz,3,1,IF] * uCol[I,J,iz+1] + 
      dXdxI[I,J,2,Iz,3,2,IF] * vCol[I,J,iz+1] + dXdxI[I,J,2,Iz,3,3,IF] * wCol[I,J,iz+1]) +
      RhoCol[I,J,iz+2] * (dXdxI[I,J,1,Iz+1,3,1,IF] * uCol[I,J,iz+2] + 
      dXdxI[I,J,2,Iz+1,3,2,IF] * vCol[I,J,iz+2] + dXdxI[I,J,1,Iz+1,3,3,IF] * wCol[I,J,iz+1])
    cL = cCol[I,J,iz+1]
    cR = cCol[I,J,iz+2]
    Flux = 0.25 * ((abs(wCon) + wCon) * cL + (-abs(wCon) + wCon) * cR)
    @inbounds @atomic FCol[I,J,iz+1] += -Flux
    @inbounds @atomic FCol[I,J,iz+2] += Flux
  end 

  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds tempx = -RhoCol[I,J,iz+1] * cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,1,1,IF] + dXdxI[I,J,2,Iz,1,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,1,2,IF] + dXdxI[I,J,2,Iz,1,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,1,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,1,3,IF] * wCol[I,J,iz+1])
    @inbounds tempy = -RhoCol[I,J,iz+1] * cCol[I,J,iz+1] * ((dXdxI[I,J,1,Iz,2,1,IF] + dXdxI[I,J,2,Iz,2,1,IF]) * uCol[I,J,iz+1] +
      (dXdxI[I,J,1,Iz,2,2,IF] + dXdxI[I,J,2,Iz,2,2,IF]) * vCol[I,J,iz+1] +
      dXdxI[I,J,1,Iz,2,3,IF] * wCol[I,J,iz] + dXdxI[I,J,2,Iz,2,3,IF] * wCol[I,J,iz+1])
    for k = 1 : N
      @inbounds @atomic FCol[k,J,iz+1] += D[k,I] * tempx
      @inbounds @atomic FCol[I,k,iz+1] += D[k,J] * tempy
    end
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz 
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind] += FCol[I,J,iz+1] / M[Iz,ind]
    if iz == 1 && Iz >  1
      @inbounds @atomic F[Iz-1,ind] += FCol[I,J,iz] / M[Iz-1,ind]
    end
    if iz == ColumnTilesDim && Iz <  Nz
      @inbounds @atomic F[Iz+1,ind] += FCol[I,J,iz+2] / M[Iz+1,ind]
    end
  end
end

@kernel function DivRhoTrUpwind3Kernel!(F,@Const(U),@Const(D),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob))

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vConCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds cCol[I,J,iz] = U[Iz,ind,5] / U[Iz,ind,1]
    @views @inbounds (uCon, vCon) = Contra12(-U[Iz,ind,1],U[Iz,ind,2],U[Iz,ind,3],dXdxI[1:2,1:2,:,ID,Iz,IF])
    @inbounds uConCol[I,J,iz] = uCon
    @inbounds vConCol[I,J,iz] = vCon
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz < Nz 
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds ind = Glob[ID,IF]
    @inbounds cL = cCol[I,J,iz]
    if iz < ColumnTilesDim 
      @inbounds cR = cCol[I,J,iz+1]
    else
      Izp1 = min(Iz + 1, Nz)
      @inbounds cR = U[Izp1,ind,5] / U[Izp1,ind,1]
    end  
        
    if iz > 1
      @inbounds cLL = cCol[I,J,iz-1]
    else
      Izm1 = max(Iz - 1,1)
      @inbounds cLL = U[Izm1,ind,5] / U[Izm1,ind,1]
    end
    if iz < ColumnTilesDim - 1
      @inbounds cRR = cCol[I,J,iz+2]
    else
      Izp2 = min(Iz + 2, Nz)
      @inbounds cRR = U[Izp2,ind,5] / U[Izp2,ind,1]
    end

    @views @inbounds wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    @inbounds JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    @inbounds JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    @inbounds JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    @inbounds JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = 0.25 * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @inbounds @atomic F[Iz,ind,5] += -Flux / M[Iz,ind]
    @inbounds @atomic F[Iz+1,ind,5] += Flux / M[Iz+1,ind]
    Flux = 0.5 * wCon
    @inbounds @atomic F[Iz,ind,1] += -Flux / M[Iz,ind]
    @inbounds @atomic F[Iz+1,ind,1] += Flux / M[Iz+1,ind]
  end 

  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds DivRho = D[I,1] * uConCol[1,J,iz] 
    @inbounds DivRho += D[J,1] * vConCol[I,1,iz] 
    @inbounds DivRhoTr = D[I,1] * uConCol[1,J,iz] * cCol[1,J,iz] 
    @inbounds DivRhoTr += D[J,1] * vConCol[I,1,iz] * cCol[I,1,iz]
    for k = 2 : N
      @inbounds DivRho += D[I,k] * uConCol[k,J,iz] 
      @inbounds DivRho += D[J,k] * vConCol[I,k,iz] 
      @inbounds DivRhoTr += D[I,k] * uConCol[k,J,iz] * cCol[k,J,iz] 
      @inbounds DivRhoTr += D[J,k] * vConCol[I,k,iz] * cCol[I,k,iz]
    end
    ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind,1] += DivRho / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,5] += DivRhoTr / M[Iz,ind]
  end
end

@kernel function DivRhoTrUpwind3Kernel!(F,@Const(U),@Const(Cache),@Const(D),@Const(DW),@Const(dXdxI),
  @Const(JJ),@Const(M),@Const(Glob),Koeff)

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF

  cCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  CacheCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  RhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  uCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  vCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  wCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FTrCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  FRhoCol = @localmem eltype(F) (N,N, ColumnTilesDim)
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds CacheCol[I,J,iz] = Cache[Iz,ind]
    @inbounds wCol[I,J,iz] = U[Iz,ind,4]
    @inbounds RhoCol[I,J,iz] = U[Iz,ind,1]
    @inbounds cCol[I,J,iz] = U[Iz,ind,5] / RhoCol[I,J,iz]
    @inbounds uCol[I,J,iz] = U[Iz,ind,2]
    @inbounds vCol[I,J,iz] = U[Iz,ind,3]
    @inbounds FRhoCol[I,J,iz] = 0
    @inbounds FTrCol[I,J,iz] = 0
  end
  @synchronize
  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz < Nz 
    ID = I + (J - 1) * N  
    @inbounds ind = Glob[ID,IF]
    @inbounds ind = Glob[ID,IF]
    @inbounds cL = cCol[I,J,iz]
    @inbounds cR = cCol[I,J,iz+1]
    if iz > 1
      @inbounds cLL = cCol[I,J,iz-1]
    else
      Izm1 = max(Iz - 1,1)
      @inbounds cLL = U[Izm1,ind,5] / U[Izm1,ind,1]
    end
    if iz < ColumnTilesDim - 1
      @inbounds cRR = cCol[I,J,iz+2]
    else
      Izp2 = min(Iz + 2, Nz)
      @inbounds cRR = U[Izp2,ind,5] / U[Izp2,ind,1]
    end

    @views @inbounds wCon = Contra3(U[Iz:Iz+1,ind,1],U[Iz:Iz+1,ind,2],U[Iz:Iz+1,ind,3],
      U[Iz,ind,4],dXdxI[3,:,:,ID,Iz:Iz+1,IF])

    Izm1 = max(Iz - 1,1)
    Izp2 = min(Iz + 2, Nz)
    @inbounds JLL = JJ[ID,1,Izm1,IF] + JJ[ID,2,Izm1,IF]
    @inbounds JL = JJ[ID,1,Iz,IF] + JJ[ID,2,Iz,IF]
    @inbounds JR = JJ[ID,1,Iz+1,IF] + JJ[ID,2,Iz+1,IF]
    @inbounds JRR = JJ[ID,1,Izp2,IF] + JJ[ID,2,Izp2,IF]
    cFL, cFR = RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR) 
    Flux = 0.25 * ((abs(wCon) + wCon) * cFL + (-abs(wCon) + wCon) * cFR)
    @inbounds @atomic F[Iz,ind,5] += -Flux / M[Iz,ind]
    @inbounds @atomic F[Iz+1,ind,5] += Flux / M[Iz+1,ind]
    Flux = 0.5 * wCon
    @inbounds @atomic F[Iz,ind,1] += -Flux / M[Iz,ind]
    @inbounds @atomic F[Iz+1,ind,1] += Flux / M[Iz+1,ind]
  end 

  if Iz <= Nz
    ID = I + (J - 1) * N  
    Dxc = 0
    Dyc = 0
    for k = 1 : N
      @inbounds Dxc = Dxc + D[I,k] * CacheCol[k,J,iz]
      @inbounds Dyc = Dyc + D[J,k] * CacheCol[I,k,iz]
    end
    
    @views @inbounds (GradDx, GradDy) = Grad12(RhoCol[I,J,iz],Dxc,Dyc,dXdxI[1:2,1:2,:,ID,Iz,IF],JJ[ID,:,Iz,IF])
    @views @inbounds (tempx, tempy) = Contra12(-Koeff,GradDx,GradDy,dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @inbounds @atomic FTrCol[k,J,iz] += DW[k,I] * tempx
      @inbounds @atomic FTrCol[I,k,iz] += DW[k,J] * tempy
    end

    @views @inbounds (tempxRho, tempyRho) = Contra12(-RhoCol[I,J,iz],uCol[I,J,iz],vCol[I,J,iz],dXdxI[1:2,1:2,:,ID,Iz,IF])
    for k = 1 : N
      @inbounds @atomic FRhoCol[k,J,iz] += D[k,I] * tempxRho
      @inbounds @atomic FRhoCol[I,k,iz] += D[k,J] * tempyRho
    end
    @inbounds tempxTr = tempxRho * cCol[I,J,iz]
    @inbounds tempyTr = tempyRho * cCol[I,J,iz]
    for k = 1 : N
      @inbounds @atomic FTrCol[k,J,iz] += D[k,I] * tempxTr
      @inbounds @atomic FTrCol[I,k,iz] += D[k,J] * tempyTr
    end
  end
  @synchronize

  Iz = (gz - 1) * ColumnTilesDim + iz
  if Iz <= Nz 
    ID = I + (J - 1) * N  
    ind = Glob[ID,IF]
    @inbounds @atomic F[Iz,ind,1] += FRhoCol[I,J,iz] / M[Iz,ind]
    @inbounds @atomic F[Iz,ind,5] += FTrCol[I,J,iz] / M[Iz,ind]
  end
end

# we may be hitting a slow path:
# https://stackoverflow.com/questions/14687665/very-slow-stdpow-for-bases-very-close-to-1
fast_powGPU(x::FT, y::FT) where {FT <: AbstractFloat} = exp(y * log(x))


@inline function Derivative(I,J,N,c,D)
  @inbounds DXc = D[I,1] * c[1,J]
  @inbounds DYc = D[J,1] * c[I,1]
  for k = 2 : N
    @inbounds DXc += D[I,k] * c[k,J]
    @inbounds DYc += D[J,k] * c[I,k]
  end
  return DXc, DYc
end  

@inline function PressureGPU(RhoTh,Phys)
  Phys.p0 * fast_powGPU(Phys.Rd * RhoTh / Phys.p0, 1 / (1 - Phys.kappa))
end  

@inline function Contra12(Rho,u,v,dXdxI)
  @inbounds uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * v)
  @inbounds vCon = Rho * ((dXdxI[2,1,1] + dXdxI[2,1,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v)
  return uCon, vCon
end

@inline function Contra12(u,v,dXdxI)
  @inbounds uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * v
  @inbounds vCon = (dXdxI[2,1,1] + dXdxI[2,1,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  return uCon, vCon
end

@inline function Grad12(Rho,u,v,dXdxI,J)
  @inbounds uCon = Rho * ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v) / (J[1] + J[2])
  @inbounds vCon = Rho * ((dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v) / (J[1] + J[2])
  return uCon, vCon
end

@inline function Grad12(u,v,dXdxI,J)
  @inbounds uCon = ((dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v) / (J[1] + J[2])
  @inbounds vCon = ((dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v) / (J[1] + J[2])
  return uCon, vCon
end

@inline function Grad12(u,v,dXdxI)
  @inbounds uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * u +
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v
  @inbounds vCon = (dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v
  return uCon, vCon
end

@inline function Grad3(u,v,dXdxI)
  @inbounds wCon1 = dXdxI[1,3,1] * u + dXdxI[2,3,1] * v
  @inbounds wCon2 = dXdxI[1,3,2] * u + dXdxI[2,3,2] * v
  return wCon1, wCon2
end

@inline function Curl12(u,v,dXdxI)
  @inbounds uCon = (dXdxI[1,1,1] + dXdxI[1,1,2]) * v -
  (dXdxI[1,2,1] + dXdxI[1,2,2]) * u 
  @inbounds vCon = (dXdxI[2,1,1] + dXdxI[2,1,2]) * v -
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * u 
  return uCon, vCon
end

@inline function Rot12(u,v,dXdxI)
  @inbounds uCon = (dXdxI[1,2,1] + dXdxI[1,2,2]) * u +
  (dXdxI[2,2,1] + dXdxI[2,2,2]) * v 
  @inbounds vCon = -(dXdxI[1,1,1] + dXdxI[1,1,2]) * u -
  (dXdxI[2,1,1] + dXdxI[2,1,2]) * v 
  return uCon, vCon
end
 
@inline function Contra3(Rho,u,v,w,dXdxI)
  wCon = Rho[1] * (dXdxI[1,2,1] * u[1] + dXdxI[2,2,1] * v[1] + dXdxI[3,2,1] * w) + 
    Rho[2] * (dXdxI[1,1,2] * u[2] + dXdxI[2,1,2] * v[2] + dXdxI[3,1,2] * w)
end
  
@inline function RecU4(cLL,cL,cR,cRR,JLL,JL,JR,JRR)

  kR = (JL / (JL + JR)) * ((JLL + JL) / (JLL + JL + JR))
  kL = -(JL / (JLL + JL)) * (JR / (JLL + JL + JR))
  cFL = kL * cLL + (1 - kL - kR)*cL + kR * cR

  kL = (JR / (JR + JL)) * ((JRR + JR)/(JL + JR + JRR))
  kR = -(JR /(JRR + JR)) *(JL /(JL + JR + JRR))
  cFR = kL * cL + (1 - kL - kR) * cR + kR * cRR
 
  return (cFL,cFR)
end

@kernel function uvwFunCKernel!(Profile,u,v,w,time,@Const(Glob),@Const(X),Param,Phys)

  gi, gz, gF = @index(Group, NTuple)
  I, iz   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = 0.5 * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = 0.5 * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = 0.5 * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,uP,vP,_ = Profile(xS,time)
    @inbounds u[Iz,ind] = uP
    @inbounds v[Iz,ind] = vP
  end
  if Iz <= Nz - 1
    ind = Glob[I,IF]
    x1 = 0.5 * (X[I,2,1,Iz,IF] + X[I,1,1,Iz+1,IF])
    x2 = 0.5 * (X[I,2,2,Iz,IF] + X[I,1,2,Iz+1,IF])
    x3 = 0.5 * (X[I,2,3,Iz,IF] + X[I,1,3,Iz+1,IF])
    xS = SVector{3}(x1, x2 ,x3)
    @inbounds _,_,_,w[Iz,ind] = Profile(xS,time)
  end
end

@kernel function RhouvFunCKernel!(Profile,Rho,u,v,time,@Const(Glob),@Const(X),Param,Phys)

  gi, gz, gF = @index(Group, NTuple)
  I, iz   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz
    ind = Glob[ID,IF]
    x1 = 0.5 * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = 0.5 * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = 0.5 * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,uP,vP,_ = Profile(xS,time)
    @inbounds Rho[Iz,ind] = RhoP
    @inbounds u[Iz,ind] = uP
    @inbounds v[Iz,ind] = vP
  end
end

@kernel function RhoFunCKernel!(Profile,Rho,time,@Const(Glob),@Const(X),Param,Phys)

  gi, gz, gF = @index(Group, NTuple)
  I, iz   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = 0.5 * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = 0.5 * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = 0.5 * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ = Profile(xS,time)
    @inbounds Rho[Iz,ind] = RhoP
  end
end

@kernel function TrFunCKernel!(Profile,Tr,time,@Const(Glob),@Const(X),Param,Phys)

  gi, gz, gF = @index(Group, NTuple)
  I, iz   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = 0.5 * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = 0.5 * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = 0.5 * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,TrP = Profile(xS,time)
    @inbounds Tr[Iz,ind] = RhoP * TrP
  end
end

@kernel function ThFunCKernel!(Profile,Th,time,@Const(Glob),@Const(X),Param,Phys)

  gi, gz, gF = @index(Group, NTuple)
  I, iz   = @index(Local, NTuple)
  _,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[2]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[2]
  NF = @uniform @ndrange()[3]

  @uniform ColumnTiles = (div(Nz - 1, ColumnTilesDim) + 1) * NF
  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz
    ind = Glob[I,IF]
    x1 = 0.5 * (X[I,1,1,Iz,IF] + X[I,2,1,Iz,IF])
    x2 = 0.5 * (X[I,1,2,Iz,IF] + X[I,2,2,Iz,IF])
    x3 = 0.5 * (X[I,1,3,Iz,IF] + X[I,2,3,Iz,IF])
    xS = SVector{3}(x1, x2 ,x3)
    RhoP,_,_,_ ,ThP = Profile(xS,time)
    @inbounds Th[Iz,ind] = RhoP * ThP
  end
end


@kernel function ComputeFunFKernel!(Profile,w,time,@Const(Glob),@Const(X),Param)

  gi, gj, gz, gF = @index(Group, NTuple)
  I, J, iz   = @index(Local, NTuple)
  _,_,_,IF = @index(Global, NTuple)

  ColumnTilesDim = @uniform @groupsize()[3]
  N = @uniform @groupsize()[1]
  Nz = @uniform @ndrange()[3]
  NF = @uniform @ndrange()[4]

  Iz = (gz - 1) * ColumnTilesDim + iz

  if Iz <= Nz - 1 
    ind = Glob[ID,IF]
    x1 = 0.5 * (X[I,J,2,1,Iz,IF] + X[I,J,1,1,Iz+1,IF])
    x2 = 0.5 * (X[I,J,2,2,Iz,IF] + X[I,J,1,2,Iz+1,IF])
    x3 = 0.5 * (X[I,J,2,3,Iz,IF] + X[I,J,1,3,Iz+1,IF])
    xS = SVector{3}(x1, x2 ,x3)
    _,_,_,wP = Profile(xS,time)
    @inbounds w[Iz,ind] = wP
  end
end

function FcnAdvectionGPU!(F,U,time,FE,Metric,Phys,Cache,Global,Param,Profile)

  backend = get_backend(F)
  FT = eltype(F)
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  N = FE.OrdPoly+1
  Nz = size(F,1)
  NF = size(Glob,3)
  Koeff = Global.Model.HyperDDiv
  Temp1 = Cache.Temp1


# State vector
  @views Rho = U[:,:,1]
  @views u = U[:,:,2]
  @views v = U[:,:,3]
  @views w = U[:,:,4]
# Cache
  @views CacheF = Temp1[:,:,1:5]
# Ranges
  NzG = min(div(512,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KuvwFunCKernel! = uvwFunCKernel!(backend, group)
  KDivRhoGradKernel! = DivRhoGradKernel!(backend, group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, group)

  KuvwFunCKernel!(Profile,u,v,w,time,Glob,X,Param,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  @. CacheF = 0
  KHyperViscKernel!(CacheF,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
  @. F = 0
  KDivRhoTrUpwind3Kernel!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,Koeff,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

@kernel function PressureKernel!(p,@Const(U),Phys)
  jz, jG = @index(Group, NTuple)
  Iz,IC = @index(Global, NTuple)

  NumG = @uniform @ndrange()[2]

  if IC <= NumG
    p[Iz,IC] = PressureGPU(U[Iz,IC,5],Phys)  
  end  
end

function FcnPrepareGPU!(U,FE,Metric,Phys,Cache,Global,Param,DiscType)
  backend = get_backend(U)
  FT = eltype(U)

  Nz = size(U,1)
  NumG = size(U,2)

  NG = div(512,Nz)
  group = (Nz, NG)
  ndrange = (Nz, NumG)
  @views p = Cache.AuxG[:,:,1]

  KPressureKernel! = PressureKernel!(backend,group)

  KPressureKernel!(p,U,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)
end

function FcnGPU!(F,U,FE,Metric,Phys,Cache,Global,Param,DiscType)

  backend = get_backend(F)
  FT = eltype(F)
  @. F = 0
  Glob = FE.Glob
  DS = FE.DS
  DW = FE.DW
  M = FE.M
  dXdxI = Metric.dXdxI
  X = Metric.X
  J = Metric.J
  DoF = FE.DoF
  N = size(FE.DS,1)
  Nz = size(F,1)
  NF = size(Glob,2)
  Koeff = Global.Model.HyperDDiv
  Temp1 = Cache.Temp1

  KoeffCurl = Global.Model.HyperDCurl
  KoeffGrad = Global.Model.HyperDGrad
  KoeffDiv = Global.Model.HyperDDiv


# State vector
  @views Rho = U[:,:,1]
  @views u = U[:,:,2]
  @views v = U[:,:,3]
  @views w = U[:,:,4]
  @views RhoTr = U[:,:,5]
# Cache
  @views CacheF = Temp1[:,:,1:6]
  @views FRho = F[:,:,1]
  @views FRhoTr = F[:,:,5]
  @views p = Cache.AuxG[:,:,1]
# Ranges
  NzG = min(div(512,N*N),Nz)
  group = (N, N, NzG, 1)
  ndrange = (N, N, Nz, NF)

  KRhoGradKinKernel! = RhoGradKinKernel!(backend,group)
  KGradKernel! = GradKernel!(backend,group)
  KDivRhoGradKernel! = DivRhoGradKernel!(backend, group)
  KHyperViscKernel! = HyperViscKernel!(backend, group)
  KHyperViscKernelKoeff! = HyperViscKernelKoeff!(backend, group)
  KDivRhoTrUpwind3Kernel! = DivRhoTrUpwind3Kernel!(backend, group)
  KMomentumKernel! = MomentumKernel!(backend, group)

  @. CacheF = 0
  @views MRho = CacheF[:,:,6]
  KHyperViscKernel!(CacheF,MRho,U,DS,DW,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  @. F = 0
  KHyperViscKernelKoeff!(F,U,CacheF,DS,DW,dXdxI,J,M,Glob,KoeffCurl,KoeffGrad,KoeffDiv,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  KGradKernel!(F,U,p,DS,dXdxI,J,M,MRho,Glob,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  KMomentumKernel!(F,U,DS,dXdxI,MRho,M,Glob,Phys,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

  KDivRhoTrUpwind3Kernel!(F,U,DS,dXdxI,J,M,Glob,ndrange=ndrange)
  KernelAbstractions.synchronize(backend)

end

