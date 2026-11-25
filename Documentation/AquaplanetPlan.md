# Plan: aquaplanet simulations with CGSEM

November 2026. Builds on the boundary-layer physics described in
[CGSEMTurbulenceSurface.md](CGSEMTurbulenceSurface.md). No code has been written for this plan
yet.

## 1 Goal

Multi-year aquaplanet climates with the CGSEM dynamical core on GPUs, in a hierarchy of
increasing physical complexity, each step compared with a published reference. The
boundary-layer physics (surface fluxes, TKE) is tested first in a configuration where all other
physics is fixed, so that differences can be attributed.

## 2 What exists

| Component | State |
|---|---|
| Dynamical core | CGSEM on the cubed sphere, moist (ρ, u, v, w, ρθ, ρq_v, ρq_c), Float32 and Float64, Rosenbrock/W with vertically implicit acoustics, about 20 SYPD dry baroclinic wave on one A100 (32 panels, OrdPoly 3, 64 levels, Δt = 150 s) |
| Surface | MOST surface layer, Charnock roughness over sea, drag, sensible and latent heat fluxes, `HeldSuarezMoistSurface` SST |
| Mixing | implicit vertical diffusion with a simple K profile, constant K, or the TKE closure |
| Moisture | `SimpleMicrophysicsPot`: saturation adjustment with immediate fall-out of a fraction of the condensate (large-scale condensation); `OneMomentMicrophysicsMoistEquil` with sedimentation |
| Forcing | `HeldSuarezMoistExample`: Held–Suarez Newtonian relaxation of T and Rayleigh friction below σ = 0.7 |
| Jobs | `Jobs/NHSphere/CG/HeldSuarezMoistSphere_32Elem` (10 days), `FriersonSphere_32Elem` (dry Held–Suarez parameters, not yet the Frierson physics) |
| Averaging | `--StartAverageDays` in the driver |

Missing for aquaplanet climates: radiation, a surface energy budget (slab ocean), a convection
scheme (or enough resolution to do without), moist buoyancy in the TKE closure, budget
diagnostics.

## 3 Hierarchy

### Stage A: moist Held–Suarez aquaplanet (Thatcher and Jablonowski 2016)

The SST of `HeldSuarezMoistSurface` (T_s = 29 K · exp(−φ²/(2 · 26°²)) + 271 K) is already the
TJ16 SST. TJ16 keeps the Newtonian relaxation of temperature, replaces the Rayleigh friction of
Held–Suarez by bulk surface fluxes and the boundary-layer diffusion of Reed and Jablonowski
(2012), and uses large-scale condensation only.

Work:

1. (done) The Rayleigh friction of `HeldSuarezMoistExample` is only applied without surface
   drag (`--SurfaceFluxMom=false`).
2. (done) `--SurfaceScheme="Bulk"`: Reed–Jablonowski coefficients C_d = 7·10⁻⁴ + 6.5·10⁻⁵ |v|
   (2·10⁻³ above 20 m/s), C_E = C_H = 1.1·10⁻³, and their K profile K = C_E |v₁| z₁ below
   850 hPa (`BulkKoefficient`; one K for momentum and scalars, RJ12 use C_d for momentum).
   Job: `Jobs/NHSphere/CG/TJ16Sphere_32Elem`.
3. Run TJ16 with these settings for 2–3 years, average after the first 6 months.
4. Compare with TJ16: zonal-mean u, T, q, relative humidity, precipitation; global mean
   precipitation and evaporation.
5. Replace bulk + K profile by MOST + TKE (one change at a time) and document the differences.

Why first: radiation and SST are fixed, so the run isolates surface fluxes, mixing and
condensation, and a published multi-model reference exists.

### Stage B: gray-radiation aquaplanet (Frierson et al. 2006)

Replaces the Newtonian relaxation with a two-stream gray radiation scheme (prescribed
longwave optical depth, shortwave absorbed at the surface) and a slab mixed-layer ocean
(heat capacity of a few metres of water); large-scale condensation only, later the simplified
Betts–Miller convection of Frierson (2007).

Work:

1. Gray radiation as a column kernel (one thread per column, upward and downward sweeps), a few
   hundred lines; heating rates into ρθ.
2. Slab ocean: prognostic T_s per column from net surface radiation, sensible and latent heat.
3. Rename the current `FriersonSphere` problem or give it the real Frierson setup.
4. Compare with Frierson et al. (2006): Hadley cell, jet position, precipitation and
   extratropical eddy kinetic energy.
5. Optional: simplified Betts–Miller convection.

Why second: the first fully interactive climate (radiation–surface–circulation) at low cost
and with a clean reference.

### Stage C: APE aquaplanet with full physics (Neale and Hoskins 2000)

Fixed QOBS SST (27 °C · (1 − sin²(3φ/2)) for |φ| < 60°, 0 °C poleward), perpetual equinox,
full radiation, clouds and convection; reference: Aqua-Planet Experiment (Blackburn et al.
2013).

Work:

1. Radiation: RRTMGP.jl (CliMA, runs on GPUs) through a thin interface; cloud optics from the
   microphysics.
2. Convection: a mass-flux scheme, or (at ≲ 10 km) none.
3. Cloud fraction: diagnostic from relative humidity, or binary from condensate.

Stage C is a separate project; the decision depends on the target resolution.

### Stage D (alternative to C): storm-resolving aquaplanet

At a few kilometres no convection scheme is needed: one-moment microphysics, TKE boundary
layer, RRTMGP. Cost estimate from Section 5 (about 10 SYPD at 100 km on one A100): halving
the grid spacing costs about a factor 8 (two horizontal directions and the time step), so
25 km is about 0.15 SYPD and 5 km about 0.001 SYPD per GPU. A storm-resolving run therefore
needs many GPUs and is limited to months.

## 4 Boundary-layer work needed before Stage A step 5

From the limitations in [CGSEMTurbulenceSurface.md](CGSEMTurbulenceSurface.md), Section 6:

| Item | Why it matters for aquaplanets |
|---|---|
| 1/Pr in CT (done) | 26 % too small latent and sensible fluxes bias evaporation and precipitation directly |
| θ₁ − θ_s instead of T₁ − T_s (done) | for z₁ ≈ 25–50 m a bias of 0.25–0.5 K, comparable to the air–sea difference over warm SST |
| Moist buoyancy (done for unsaturated air) | N² and the surface buoyancy flux with θ_v (and the saturated N² in clouds); otherwise the marine boundary layer is too stable |
| Implicit surface heat and moisture fluxes | Δt = 150–300 s with Δz₁ ≈ 30 m: CT·U_a·Δt/Δz₁ approaches 0.1–0.2, more over strong winds |
| Gustiness (done, constant `--Gustiness`) | a minimum wind of about 1 m/s (Beljaars 1995) instead of 0.1 m/s keeps fluxes in the calm subtropics |
| Mixing length | λ of 100 m is set for the rotor case; for the global boundary layer scale it with the boundary-layer depth or use λ ≈ 150 m |
| Water and energy budgets | global E − P and the energy budget as running diagnostics; needed to trust multi-year means |

## 5 Configuration

| | Stage A/B |
|---|---|
| Grid | cubed sphere, 32 panels, OrdPoly 3 (≈ 1°, 100 km) |
| Vertical | 64 levels, top 30 km, Exp stretching with Δz₁ ≈ 30–50 m (10–15 levels below 2 km) |
| Time step | 150 s, Rosenbrock SSP-Knoth |
| Precision | Float32 for production, Float64 for the first checks |
| Tracers | q_v, q_c (`SimpleMicrophysicsPot`), TKE |
| Length | 3 years, average years 2–3 |
| Cost | about 10 SYPD on one A100 (physics overhead assumed ≤ 2×), i.e. a 3-year run in a few hours |

## 6 Diagnostics

- Zonal and time means: u, v, ω, T, θ, q_v, RH, precipitation, evaporation, sensible heat flux,
  TKE, K, boundary-layer height.
- Global means: precipitation = evaporation in equilibrium, energy budget at the top (Stage B).
- Mass stream function (Hadley cell), eddy kinetic energy, ITCZ position and width.
- Output as zonal means on the GPU node, plotted with Python batch jobs on Levante (no data
  transfer to the Mac).

## 7 Order of work

1. Boundary-layer fixes of Section 4 (Pr, θ difference, gustiness, moist buoyancy), column
   tests as in [CGSEMTurbulenceSurface.md](CGSEMTurbulenceSurface.md), Section 5.1.
2. Zonal-mean and budget diagnostics.
3. Stage A with the TJ16 settings; agreement with TJ16 is the acceptance test.
4. Stage A with MOST + TKE.
5. Implicit surface scalar fluxes, if the time step requires them.
6. Stage B (gray radiation, slab ocean).
7. Decision on Stage C or D.

## 8 Decisions for you

1. Start with TJ16 (Stage A) as proposed, or go directly to Frierson (Stage B)?
2. Fixed SST or slab ocean in Stage B?
3. Radiation beyond gray: RRTMGP.jl, or something else?
4. Convection: large-scale condensation only, simplified Betts–Miller, or high resolution
   without a scheme?
5. Should DGSEM follow each stage, or only after Stage B?

## References

- Beljaars, A. C. M., 1995: The parametrization of surface fluxes in large-scale models under
  free convection. Quart. J. Roy. Meteor. Soc., 121, 255–270.
- Blackburn, M., et al., 2013: The Aqua-Planet Experiment (APE): control SST simulation.
  J. Meteor. Soc. Japan, 91A, 17–56.
- Frierson, D. M. W., I. M. Held, and P. Zurita-Gotor, 2006: A gray-radiation aquaplanet moist
  GCM. Part I. J. Atmos. Sci., 63, 2548–2566.
- Frierson, D. M. W., 2007: The dynamics of idealized convection schemes and their effect on
  the zonally averaged tropical circulation. J. Atmos. Sci., 64, 1959–1976.
- Neale, R. B., and B. J. Hoskins, 2000: A standard test for AGCMs including their physical
  parametrizations: I: The proposal. Atmos. Sci. Lett., 1, 101–107.
- Reed, K. A., and C. Jablonowski, 2012: Idealized tropical cyclone simulations of intermediate
  complexity. J. Adv. Model. Earth Syst., 4, M04001.
- Thatcher, D. R., and C. Jablonowski, 2016: A moist aquaplanet variant of the Held–Suarez test
  for atmospheric model dynamical cores. Geosci. Model Dev., 9, 1263–1292.
