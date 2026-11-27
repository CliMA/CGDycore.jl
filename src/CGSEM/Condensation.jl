# Large-scale condensation as a saturation adjustment after every step
# (--TypeMicrophysics="Adjustment"), as in the moist Held-Suarez test of Thatcher and Jablonowski
# (2016): supersaturated vapour condenses and falls out at once. Condensation at constant volume
# with the internal energies of Thermodynamics (heat e_v - e_l = L - Rv*T per condensed mass), the
# condensate leaves at the temperature of the air, rho theta follows from the moist pressure law.
# Newton steps until rho_vs(T) - rho_v = (1 - f) * (rho_vs - rho_v) of the start, f = min(1, dtau/TauR):
# TauR <= dtau is the full adjustment, a larger TauR a relaxation with this time scale.
# Precip holds the column precipitation rate [kg/m^2/s] of the last step.
struct Condensation{FT,AT1,PT}
  Precip::AT1
  Phys::PT
  dtau::FT
  Frac::FT
end

function Condensation(backend,FT,NumG,Phys,dtau;TauR=0)
  Precip = KernelAbstractions.zeros(backend,FT,NumG)
  Frac = TauR > dtau ? FT(dtau / TauR) : FT(1)
  return Condensation{FT,typeof(Precip),typeof(Phys)}(Precip,Phys,FT(dtau),Frac)
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
  # heating at constant volume: dT = (L - Rv*T) / Cvml per condensed mass
  Cvml = Phys.Cvd * RhoD + Phys.Cvv * RhoV + Phys.Cpl * RhoC
  Q = L - Phys.Rv * T
  dadd = FT(1) + RhoVS * (L / (Phys.Rv * T * T) - FT(1) / T) * Q / Cvml
  # Rd*rho*theta = Rm*T*(p0/p)^kappaM with p = Rm*T, the vapour leaves Rm and Cpml:
  # dln(rho theta) = (1 - kappaM)*(dT/T - Rv/Rm) - ln(p/p0)*dkappaM, (1 - kappaM)/Cvml = 1/Cpml
  lpp0 = log(p / Phys.p0)
  dRhoTh = RhoTh * (Q / (Cpml * T) - (FT(1) - kappaM) * Phys.Rv / Rm +
    lpp0 * kappaM * (Phys.Rv / Rm - Phys.Cpv / Cpml))
  return a, dadd, dRhoTh
end

@kernel inbounds = true function CondensationKernel!(U,Precip,@Const(dz),Phys,ThPos,RhoVPos,RhoCPos,dtau,Frac)
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
        # target supersaturation of the relaxation, zero for the full adjustment
        aT = (FT(1) - Frac) * a
        # Newton steps, back steps allowed as long as the condensed mass stays positive
        for Iter = 1 : 5
          d = max((aT - a) / dadd,-Cond)
          Rho -= d
          RhoV -= d
          RhoTh += dRhoTh * d
          Cond += d
          a, dadd, dRhoTh = CondensationState(Rho,RhoTh,RhoV,RhoC,Phys)
          if abs(a - aT) <= FT(1.0e-4) * RhoV
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
    Global.Model.RhoCPos,C.dtau,C.Frac,ndrange=(NumG,))
end
CondensationAdjust!(U,Metric,Global,C::Nothing) = nothing
