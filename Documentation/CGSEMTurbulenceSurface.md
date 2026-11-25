# Turbulence and surface fluxes in CGSEM

State of branch `cgsem-pbl` (November 2026). This note describes what the code does, how to
switch it on, how it was tested, and the known limitations. The plan that led to it is in
[PBLPlan.md](PBLPlan.md); the aquaplanet plan that builds on it is in
[AquaplanetPlan.md](AquaplanetPlan.md).

## 1 Overview

The boundary-layer physics of CGSEM has three parts, all evaluated column by column:

| Part | What it computes | Code |
|---|---|---|
| Surface layer | u\*, bulk coefficients CM, CT, CH, wind speed U_a from Monin–Obukhov similarity | `src/Surfaces/surface.jl`, `solver.jl`, `functions.jl` |
| Surface fluxes | drag on u, v; sensible heat and moisture fluxes into the lowest layer | `src/Surfaces/surface.jl`, `src/CGSEM/MomentumKernel.jl` |
| Vertical mixing | eddy diffusivity K (TKE closure, constant, or simple profile) and implicit vertical diffusion of u, v, θ, tracers and TKE | `src/Models/Turbulence.jl`, `src/CGSEM/DiagnosticKernel.jl`, `Jac.jl`, `Solve.jl` |

Horizontal transport of TKE uses the normal tracer advection; there is no horizontal turbulent
diffusion apart from the hyperviscosity of the dynamical core.

## 2 Surface layer

### 2.1 Surface values

`Model.SurfaceValues` sets the surface temperature T_s and the surface water vapour density
ρ_vs per column (`SurfaceData`, positions `TSurfPos`, `RhoVSurfPos`):

- `DefaultSurface`: constant `Param.TSurf`, `Param.RhoVSurf`.
- `HeldSuarezMoistSurface`: T_s = ΔT_s exp(−φ²/(2Δφ²)) + T_s,min (ΔT_s = 29 K, T_s,min = 271 K,
  Δφ = 26°), ρ_vs saturated at T_s. Used for `HeldSuarezMoistSphere`, `BaroWaveMoistSphere`,
  `FriersonSphere`.

### 2.2 Monin–Obukhov similarity (`--SurfaceScheme="MOST"`)

`MOSurfaceFlux` evaluates the similarity laws at the centre of the lowest layer,
z₁ = Δz₁/2, with the layer values of the wind speed U_a and θ (the "conventional" scheme in the
terminology of Nishizawa and Kitamura 2018). The universal functions are those of Businger et
al. (1971) (`Businger`, a = 4.7, Pr = 0.74). The bulk Richardson number
Ri_b = g z₁ (θ₁ − θ_s) / (θ_s U_a²), clipped to [−10, 0.2], is inverted for ζ = z₁/L by Newton
iteration (`MOSTIteration`), and

- CM = κ² / (ln(z₁/z₀m) − ψ_m(ζ))²
- CT = CH = κ² / [(ln(z₁/z₀m) − ψ_m(ζ)) (ln(z₁/z₀h) − ψ_h(ζ))]
- u\* = κ U_a / (ln(z₁/z₀m) − ψ_m(ζ)) = √CM · U_a

CM, CT, CH are dimensionless bulk coefficients: every flux is coefficient × U_a × difference.
This is the same as the u\*θ\* form, F = −u\* · κ/(ln(z₁/z₀h) − ψ_h) · Δθ.

U_a is the speed of the velocity tangential to the surface at the lowest level, bounded below
by 0.1 m/s (`UzCoefficientGPU`).

Roughness (`LandUseData`, two classes, `--LandClass`):

| Class | z₀m, z₀h |
|---|---|
| 1 land | `--z0M`, `--z0H` (default 0.01 m) |
| 2 sea (default) | Charnock-type from u\*: z₀ = exp((2.7u\* − 14.4)/(u\* + 1.39)), iterated with ζ (`MOSTSeaIteration`); `--z0M`/`--z0H` are ignored |

### 2.3 Surface fluxes

Momentum (`--SurfaceFluxMom`, `SurfaceFluxMomentumKernel!`):

  ∂u/∂t (layer 1) −= DragFac · CM · U_a · u / Δz₁

`DragFac` (per column, default 1) allows free-slip regions; it is used by the rotor case.

Heat and moisture (`--SurfaceFlux`, `SurfaceFlux`):

  F_θ = −CT · U_a · (T₁ − T_s),  F_v = −CH · U_a · (ρ_v1 − ρ_vs)

added to the lowest layer as F/Δz₁. In the moist version the tendency of ρθ also contains the
change of θ from the moisture flux.

## 3 Vertical mixing

### 3.1 Eddy diffusivity

`Model.Eddy` selects K. It is stored as ρK at the layer interfaces (`Cache.KV[1,iz,ic]` is the
interface iz+½; the top interface iz = Nz is zero):

| Option | K |
|---|---|
| `--Turbulence=true` | TKE closure (3.2) |
| `--KConst > 0` | constant K (tests) |
| otherwise | `SimpleKoefficient`: K = CE · U_a · Δz₁/2, decaying as exp(−((p_pbl − p)/p_strato)²) above p_pbl = 850 hPa (Reed and Jablonowski 2012) |

### 3.2 TKE closure (`TKEModel`, `--Turbulence=true --TkePos=k`)

A one-and-a-half order closure with prognostic ρe (Mellor–Yamada level 2.5 constants):

| Quantity | Formula | Default |
|---|---|---|
| eddy diffusivity | K_m = K_h = c_m l √e | c_m = 0.554 (Pr_t = 1) |
| dissipation | ε = c_ε e^{3/2} / l | c_ε = c_m³ ≈ 0.17 |
| mixing length | l = κz / (1 + κz/λ), in stable layers ≤ 0.75 √(e/N²), ≥ 0.1 m | λ = `--MixingLength` (100 m) |
| production | P = K (S² − N²) at the interfaces | |
| surface production | u\*³/(κ z₁) + g F_θ/T₁ | with `--SurfaceFluxMom` |
| minimum | e ≥ 10⁻⁴ m²/s² | |

With c_ε = c_m³ the neutral surface layer gives K = κu\*z and e = u\*²/c_m².
S² and N² use differences between neighbouring layers; N² uses the dry potential temperature.

`TkeColumnKernel!` (one thread per column) computes K at the interfaces, the source
ρ(P − ε) per layer (mean of the production at the two faces) and ∂ε/∂e for the implicit solve.

### 3.3 Vertical diffusion

Scalars (θ, tracers, TKE) are diffused as ∂(ρq)/∂t = ∂/∂z(ρK ∂q/∂z), the horizontal velocity as
∂u/∂t = (1/ρ) ∂/∂z(ρK ∂u/∂z) (`--VerticalDiffusion`, `--VerticalDiffusionMom`). The surface flux
is the lower boundary condition; the top flux is zero.

### 3.4 Implicit treatment (`--JacVerticalDiffusion=true`)

Inside the Rosenbrock/W stages the vertical diffusion is treated implicitly by approximate
factorisation: a tridiagonal diffusion solve per column comes before the acoustic Schur solve.

- `JacTransportKernel!` assembles the diffusion Jacobian row by row, separately for scalars
  (`JDiff`, off-diagonals with ρ_j) and velocities (`JDiffMom`, with ρ_i).
- The surface drag DragFac·CM·U_a/Δz₁ (`JDrag`) and the TKE sink ∂ε/∂e (`JDiss`) are on the
  diagonal.
- `SchurSolveDiffKernel!` solves all components with the Thomas algorithm.

The surface heat and moisture fluxes are explicit.

## 4 Switching it on

| Option | Meaning |
|---|---|
| `--SurfaceScheme="MOST"` | required as soon as any of the switches below is on; `"Bulk"` gives the constant-type coefficients and K profile of Reed and Jablonowski (2012) |
| `--SurfaceFlux`, `--SurfaceFluxMom` | surface heat/moisture fluxes, surface drag |
| `--VerticalDiffusion`, `--VerticalDiffusionMom` | vertical diffusion of scalars, velocity |
| `--JacVerticalDiffusion` | implicit vertical diffusion |
| `--Turbulence --TkePos=k --MixingLength=λ` | TKE closure, TKE stored in U[:,:,:,k] |
| `--KConst=K` | constant K instead |
| `--LandClass`, `--z0M`, `--z0H` | surface type and roughness |
| `--StretchType="Geom" --dzBottom --dzTop` | geometric vertical grid: Δz₁ = dzBottom, growing to at most dzTop |

Example jobs: `Jobs/NHCart/CG/ColumnDiffusion`, `ColumnTKE` (single columns),
`Jobs/NHCart/CG/RotorCart` (2D mountain rotors).
`Examples/PBLRun.jl` runs any DriverCG job and writes column profiles (`PBLMODE=column`) or
x–z slices (`PBLMODE=slice`) instead of VTK files.

## 5 Tests

### 5.1 Column tests (Mac CPU, 200 km × 200 km box, 2 × 2 elements)

| Test | Result |
|---|---|
| Gaussian bump in u, K = 10 m²/s, 1 h | implicit Δt = 60 s matches the analytic solution: peak 1.7462 against 1.7461 |
| θ diffusion | implicit 60 s against explicit 5 s: 5·10⁻⁵ K |
| Neutral, U = 10 m/s, z₀ = 0.1 m, TKE, 2 h | u\* = κU₁/ln(z₁/z₀); log law within 1–2 % up to 90 m; K ≈ 0.85–0.9 κu\*z; e ≈ u\*²/c_m²; implicit 60 s against explicit 1 s: 10⁻³ |
| Heated / cooled surface (±5 K) | unstable: TKE up to 0.9 m²/s², K ≈ 30 m²/s; stable: shallow layer, K ≈ 1 m²/s |
| Float32 against Float64 | 10⁻⁴ |

### 5.2 Mountain rotors (Doyle and Durran 2002, Section 3)

Two-layer flow, U = 25 m/s, N = 0.025 s⁻¹ below 3 km and 0.01 s⁻¹ above, Witch of Agnesi ridge
600 m high with half-width 10 km, z₀ = 0.1 m, no surface heat flux, free slip upstream of
x = −2a. Domain 259 km at 100 m node spacing, Δz from 10 m at the ground to 100 m, top 15 km
with w damping above 10 km, one A100, Float32, Δt = 0.1 s.

Periodic in x, 3 h (68 min wall time): the flow separates on the lee slope 8–10 km behind the
crest, trapped lee waves of 7.5–8 km wavelength form, and reversed flow appears beneath the
crests in the lowest 100–150 m. The strongest rotor reaches −7 m/s at x = 13–20 km between 1.5
and 2.25 h; maximum TKE about 20 m²/s². After about 2 h the flow re-entering through the
periodic boundary disturbs the upstream profile.

Open boundaries in x (`--BoundaryWE="Open"`, `src/CGSEM/LateralBoundary.jl`):

- inflow (west): the whole state keeps its initial values (zero tendency, reset after every
  step). With a free density and the one-sided strong-form derivative the column next to the
  inflow became unstable;
- outflow (east): one ghost element column east of the outflow x_out = x0 + Lx, added by the
  driver (nx + 1 elements). Before every right-hand side and after every step the ghost nodes
  get the state mirrored at x_out; the topography in the ghost elements is mirrored as well.
  After DSS the outflow node is an interior node with zero x derivative of the whole state,
  the analogue of the outflow face of a staggered finite-volume scheme, and hyperviscosity and
  upwinding act there as in the interior.

Outflow variants that failed (4 km test box, CPU, Float64):

| Outflow | Result |
|---|---|
| strong form only ("do nothing") | blew up within 7 s |
| weak boundary term p′ = p − p₀ | stable on the CPU, 2.25 h on the GPU, then the outflow boundary layer failed |
| zero horizontal pressure gradient | outflow column stagnated, NaN after 1.75 h on the GPU |
| copy or polynomial zero gradient at the outflow node | acoustic instability within 25–60 s, independent of Δt, also in uniform flow over flat ground |
| mirrored ghost element | stable (uniform flow 120 s, rotor case below) |

A reduced rotor run on the Mac (8 threads, 160 km from x = −40 km, Δx = 200 m, Δt = 0.2 s,
hyperviscosity scaled with Δx⁴ to 4.8·10⁵ m⁴/s; with 3·10⁴ it blew up in the lee after 4 min)
reached 29 min with a quiet outflow (|w| ≤ 0.01 m/s) and a downslope windstorm of 45 m/s. The
fixed inflow, however, traps upstream-propagating long gravity waves of the start-up: |w| at
x = −39 km, z = 2.6 km grew from 0.6 m/s at 9 min to 6.6 m/s at 29 min. A ghost element at the
inflow (u, w, θ prescribed, ρ mirrored) or a relaxation zone is the next step for this case.

## 6 Known limitations

1. **Pr in CT (fixed).** CT = κ²/(V_M · Pr · V_H) with the Businger heat profile
   φ_h = Pr + 4.7ζ, consistent with the Ri_b iteration; before, the factor 1/Pr was missing and
   CT was 26 % too small. Heated column (T_s = 305 K, 2 h): CT 0.0097 → 0.0125, θ at 10 m
   302.88 → 303.13 K.
2. **Temperature difference (fixed).** Ri_b and the sensible heat flux use θ₁ − θ_s with
   θ_s = T_s (p₀/p_s)^κ and p_s = p₁ + ρ₁ g Δz₁/2; before, they used θ₁ − T_s and T₁ − T_s,
   which is wrong by g z₁/c_p and, in Ri_b, by T_s((p₀/p_s)^κ − 1) (8 K at p_s = 900 hPa).
3. **Layer averages.** The similarity laws are applied at z₁ = Δz₁/2 to layer averages. For
   z₀ = 0.1 m and Δz₁ = 10 m this underestimates u\* by about 8 % (Nishizawa and Kitamura 2018).
   We keep this conventional scheme for now.
4. **Unsaturated buoyancy only.** N² uses θ_v = θ (1 + (R_v/R_d − 1) q_v − q_c) and the surface
   buoyancy flux g (w′θ′/T₁ + (R_v/R_d − 1) w′q_v′); the surface buoyancy term is only active
   with `--SurfaceFlux=true` (before, it was also added without surface heat flux). There is no
   saturated N² in clouds and no cloud-top mixing.
5. **Explicit surface heat and moisture fluxes.** Only the drag is in the Jacobian; for thin
   lowest layers and long time steps the heat flux can become stiff.
6. **Minimum wind.** U_a has a floor of 0.1 m/s; `--Gustiness=u_g` (default 0) adds a constant
   gustiness velocity, U_a = (|v|² + u_g²)^½, for the MOST scheme. A convective gustiness
   (β w\*, Beljaars 1995) is not implemented.
7. **Mixing length** λ is a fixed number, not scaled with the boundary-layer depth.
8. **DGSEM** has none of this yet; see Section 4 of [PBLPlan.md](PBLPlan.md).

## References

- Businger, J. A., et al., 1971: Flux-profile relationships in the atmospheric surface layer.
  J. Atmos. Sci., 28, 181–189.
- Doyle, J. D., and D. R. Durran, 2002: The dynamics of mountain-wave-induced rotors.
  J. Atmos. Sci., 59, 186–201.
- Mellor, G. L., and T. Yamada, 1982: Development of a turbulence closure model for geophysical
  fluid problems. Rev. Geophys., 20, 851–875.
- Nishizawa, S., and Y. Kitamura, 2018: A surface flux scheme based on the Monin-Obukhov
  similarity for finite volume models. J. Adv. Model. Earth Syst., 10, 3159–3175.
- Reed, K. A., and C. Jablonowski, 2012: Idealized tropical cyclone simulations of intermediate
  complexity. J. Adv. Model. Earth Syst., 4, M04001.
