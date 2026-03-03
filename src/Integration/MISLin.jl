mutable struct CacheMISLinStruct{FT<:AbstractFloat,
                           AT4<:AbstractArray,
                           AT5<:AbstractArray}
  Vn::AT4
  dZn::AT4
  fV::AT4
  Sdu::AT5
  Yn::AT5
  Delta::AT4
# ROS
  k::AT5
end

function Cache(backend,FT,IntMethod::MISLinMethod,FE,M,nz,NumV) 
  NumG = FE.NumG
  NumI = FE.NumI
  nStage = IntMethod.nStage
  Vn = KernelAbstractions.zeros(backend,FT,M,nz,NumG,NumV)

  dZn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  fV = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)
  Sdu = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage)
  Yn = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,nStage+1)
  Delta = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV)

  k = KernelAbstractions.zeros(backend,FT,M,nz,NumI,NumV,IntMethod.FastMethod.nStage)
  return CacheMISLinStruct{FT,
                     typeof(Vn),
                     typeof(k)}(
    Vn,
    dZn,
    fV,
    Sdu,
    Yn,
    Delta,
    k,
  )
end

#=
function y=SplitExplicitLin(y0,dt,Time,MIS,Fcn,Jac)
global Param
ns=Param.ns;
dtInt=dt;
Y=zeros(size(y0,1),MIS.nStage+1);
FY=zeros(size(y0,1),MIS.nStage);

for iStage=1:MIS.nStage+1
  Y(:,iStage)=y0;
end
[JE,JI]=Jac(Time,y0);
for iStage=1:MIS.nStage
  for jStage=1:iStage
    Y(:,iStage+1)=MIS.alpha(iStage+1,jStage)*(Y(:,jStage)-y0)+Y(:,iStage+1);
  end
  %FY(:,iStage)=Fcn(Time+MIS.c(iStage)*dtInt,Y(:,iStage))-JE*Y(:,iStage)-JI*Y(:,iStage);
  FY(:,iStage)=Fcn(Time+MIS.c(iStage)*dtInt,Y(:,iStage));
  FSlow=zeros(size(y0,1),1);
  Delta=zeros(size(y0,1),1);
  for jStage=1:iStage
    FSlow=FSlow+MIS.betaS(iStage+1,jStage).Koeff(1)*FY(:,jStage)+MIS.gammaS(iStage+1,jStage).Koeff(1)/dtInt*(Y(:,jStage)-y0);
    Delta=Delta+MIS.betaS(iStage+1,jStage).Koeff(1)*Y(:,jStage);
  end
  Y(:,iStage+1)=Y(:,iStage+1)-Delta;
  nsLoc=ceil(ns*MIS.d(iStage+1));
  dtLoc=dtInt*MIS.d(iStage+1);
  dTau=dtLoc/nsLoc;
  switch Param.FastIntMethod
    case 'Verlet'
      Y(:,iStage+1)=VerletLin(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE+JI)+Delta;
    case 'ForBack'
      Y(:,iStage+1)=ForBackLin(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE+JI)+Delta;
    case 'ExpMatrix'
      Y(:,iStage+1)=ExpMatrixLin(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE+JI)+Delta;
    case 'IMEX'
      Y(:,iStage+1)=LinIMEX_RK(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE,JI);
    case 'EXRK'
      Y(:,iStage+1)=Lin_RK(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE+JI);
    case 'EXLSRK'
      Y(:,iStage+1)=Lin_LSRK(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE+JI);
      fprintf(['%5i'],nsLoc);
    case 'LSVerlet'
      Y(:,iStage+1)=Lin_LSVerlet(Y(:,iStage+1),FSlow,nsLoc,dTau,Time,JE+JI);
      fprintf(['%5i'],nsLoc);
  end
end
fprintf('  h= %5.3e, t1 = %5.3e\n',dt,Time);
y=Y(:,MIS.nStage+1);
end
=#


function TimeIntegration!(MIS::MISLinMethod,V,dt,Fcn,Aux,Jac,FE,Metric,Phys,Cache,JCache,Exchange,
  Global,Param,Type)

  FcnSlow, FcnFast = Fcn
  dtSlow,dtFast = dt
  Fast = MIS.FastMethod

  nStage = MIS.nStage
  k = Cache.k
  Vn = Cache.Vn
  Sdu = Cache.Sdu
  dZn = Cache.dZn
  Yn = Cache.Yn
  Delta = Cache.Delta
  @views VnI = Vn[:,:,1:size(V,3),1:size(V,4)]
  Global.TimeStepper.dtauStage = dtSlow

  @views @. Yn[:,:,:,:,1] = V
  @inbounds for iStage = 1 : nStage
    @. VnI = Yn[:,:,:,:,iStage]
    @views FcnSlow(Sdu[:,:,:,:,iStage],Vn,FE,Metric,Phys,Aux,Exchange,Global,Type)
    if iStage == 1
      @views FrozenStateLin(Vn,Aux.Aux)
    end  
    @views SlowTendencyLin!(dZn, Delta, Yn, Sdu, V, MIS.d, MIS.gamma, MIS.beta, dtSlow, iStage)
    @views @.  Yn[:,:,:,:,iStage+1] = V
    @views InitialConditionMIS!(Yn[:,:,:,:,iStage+1], Yn, V, MIS.alpha, iStage)
    @views TimeStepperFast!(Fast,dtFast,MIS.d[iStage+1]*dtSlow,Yn[:,:,:,:,iStage+1],Delta,dZn,FcnFast,Jac,FE,
       Exchange,Metric,Phys,Param,Global,Cache,JCache,Aux,Vn,k,Type)
  end
  @views @. V = Yn[:,:,:,:,end]
end  

@inline function FrozenStateLin(Vn,Aux)
  @views @. Aux[:,:,:,3] = Aux[:,:,:,1] / Vn[:,:,:,5] 
  @views @. Aux[:,:,:,4] = Vn[:,:,:,5] / Vn[:,:,:,1]
end

@inline function SlowTendencyLin!(dZn, D, Yn, Sdu, u, d, gamma, beta, dtL, stage)
  @views @. dZn = (beta[stage+1,1] / d[stage+1]) * Sdu[:,:,:,:,1]
  @views @. D = (beta[stage+1,1] / d[stage+1]) * Yn[:,:,:,:,1]
  @inbounds for j in 2 : stage
    @views @. dZn += (gamma[stage+1,j] / (dtL * d[stage+1])) * (Yn[:,:,:,:,j] - u) +
      (beta[stage+1,j] / d[stage+1]) * Sdu[:,:,:,:,j]
    @views @. D += (beta[stage+1,j] / d[stage+1]) * Yn[:,:,:,:,j]
  end
end


