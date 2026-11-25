# Large-scale condensation as a saturation adjustment after every step
# (--TypeMicrophysics="Adjustment"), as in the moist Held-Suarez test of Thatcher and Jablonowski
# (2016): supersaturated vapour condenses and falls out at once. The increments of rho, rho_v and
# rho theta per condensed mass are those of Models.SimpleMicrophysicsPot with Rain = 1 (heating at
# constant volume), applied in Newton steps until rho_v = rho_vs(T); this is the limit of an
# infinite relaxation rate and has no stability limit. Precip holds the column precipitation rate
# [kg/m^2/s] of the last step.
struct Condensation{FT,AT1,PT}
  Precip::AT1
  Phys::PT
  dtau::FT
end

function Condensation(backend,FT,NumG,Phys,dtau)
  Precip = KernelAbstractions.zeros(backend,FT,NumG)
  return Condensation{FT,typeof(Precip),typeof(Phys)}(Precip,Phys,FT(dtau))
end

# supersaturation a = rho_vs(T) - rho_v, its derivative with respect to the condensed mass and
# the increment of rho theta per condensed mass
@inline function CondensationState(Rho,RhoTh,RhoV,RhoC,Phys)
  FT = eltype(Rho)
  RhoD = Rho - RhoV - RhoC
  Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoC
  Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
  kappaM = Rm / Cpml
  p = (Phys.Rd * RhoTh / Phys.p0^kappaM)^(FT(1) / (FT(1) - kappaM))
  T = p / Rm
  RhoVS = Thermodynamics.fpws(T,Phys) / (Phys.Rv * T)
  L = Thermodynamics.LatHeatV(T,Phys)
  a = RhoVS - RhoV
  # the increments heat at constant volume: dT = L / Cvml per condensed mass
  Cvml = Phys.Cvd * RhoD + Phys.Cvv * RhoV + Phys.Cpl * RhoC
  dadd = FT(1) + RhoVS * (L / (Phys.Rv * T * T) - FT(1) / T) * L / Cvml
  # SimpleMicrophysicsPot with FPh = -1, FR = 1
  lpp0 = log(p / Phys.p0)
  dRhoTh = RhoTh * ((L / (Cpml * T) + lpp0 * kappaM * (Phys.Rv / Rm + (Phys.Cpl - Phys.Cpv) / Cpml) -
    Phys.Rv / Rm) + (FT(1) / Rho - lpp0 * kappaM * (Phys.Cpl / Cpml)))
  return a, dadd, dRhoTh
end

@kernel inbounds = true function CondensationKernel!(U,Precip,@Const(dz),Phys,ThPos,RhoVPos,RhoCPos,dtau)
  IC, = @index(Global, NTuple)
  NumG = @uniform @ndrange()[1]

  if IC <= NumG
    FT = eltype(U)
    Nz = size(U,2)
    P = FT(0)
    for Iz = 1 : Nz
      Rho = U[1,Iz,IC,1]
      RhoTh = U[1,Iz,IC,ThPos]
      RhoV = U[1,Iz,IC,RhoVPos]
      RhoC = RhoCPos > 0 ? U[1,Iz,IC,RhoCPos] : FT(0)
      Cond = FT(0)
      a, dadd, dRhoTh = CondensationState(Rho,RhoTh,RhoV,RhoC,Phys)
      if a < FT(0)
        # Newton steps, back steps allowed as long as the condensed mass stays positive
        for Iter = 1 : 5
          d = max(-a / dadd,-Cond)
          Rho -= d
          RhoV -= d
          RhoTh += dRhoTh * d
          Cond += d
          a, dadd, dRhoTh = CondensationState(Rho,RhoTh,RhoV,RhoC,Phys)
          if abs(a) <= FT(1.0e-4) * RhoV
            break
          end
        end
      end
      if Cond > FT(0)
        U[1,Iz,IC,1] = Rho
        U[1,Iz,IC,ThPos] = RhoTh
        U[1,Iz,IC,RhoVPos] = RhoV
        P += Cond * dz[Iz,IC]
      end
    end
    Precip[IC] = P / dtau
  end
end

function CondensationAdjust!(U,Metric,Global,C::Condensation)
  backend = get_backend(U)
  NumG = size(U,3)
  group = (min(Global.ParallelCom.NumberThreadGPU,NumG),)
  KCondensationKernel! = CondensationKernel!(backend,group)
  KCondensationKernel!(U,C.Precip,Metric.dz,C.Phys,Global.Model.RhoThPos,Global.Model.RhoVPos,
    Global.Model.RhoCPos,C.dtau,ndrange=(NumG,))
end
CondensationAdjust!(U,Metric,Global,C::Nothing) = nothing
