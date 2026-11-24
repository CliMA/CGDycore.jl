# Surface fluxes, vertical diffusion and TKE: state and plan

The surface-layer, vertical-diffusion and TKE code of CGDycore.jl was reviewed
in CGSEM (the only dycore that uses it). This note plans how to make it
consistent and then bring it to DGSEM. No code was changed.

## 0 Status (branch `cgsem-pbl`)

Steps 1–4 of section 3 are implemented in CGSEM:

- **Surface layer:**
  - the bulk fluxes use the wind speed |v| (new slot `UaPos`), and `uStarPos` holds u\*;
  - MOST is evaluated at Δz₁/2;
  - the momentum flux is divided by Δz₁ and scaled by a per-column drag factor
    (`LandUseData.DragFac`; 1 no-slip, 0 free-slip);
  - `zetaPos` no longer shares its slot with `hBLPos`;
  - Float32-clean literals.
- **New options:**
  - `--LandClass` (1 land with `--z0M`/`--z0H`, 2 sea with the Charnock roughness, default 2);
  - `--KConst` (constant K for tests);
  - `--MixingLength`;
  - `--StretchType="Geom"` with `--dzBottom`/`--dzTop`.
- **Eddy diffusivity:** at the interfaces (`KV[1,Iz,IC]` for Iz+1/2) for every scheme; the
  momentum kernel reads it correctly.
- **TKE:**
  - `Models.TKEModel`, a Mellor–Yamada 2.5-type closure:
    - K = c_m l √e with c_m = 0.554, and K_h = K_m;
    - dissipation c_ε e^{3/2}/l with c_ε = c_m³ = 0.17;
    - Blackadar mixing length, limited by 0.75√(e/N²);
    - production K(S² − N²), plus the surface production u\*³/(κz₁) + g w'T'/T;
  - `TkeColumnKernel!` computes it with one thread per column.
- **Implicit column solve:**
  - row-wise Jacobians: `JDiff` for the ρq scalars and `JDiffMom` for the velocities
    (the old one had swapped and shifted entries);
  - the linearised drag (`JDrag`) and the TKE sinks (`JDiss`) on the diagonal;
  - a Thomas solve per column (`SchurSolveDiffKernel!`).
- **Driver fixes:**
  - `AgnesiCartX`/`AgnesiHill` now use P1–P3 (before, `AgnesiCartX` gave flat ground);
  - the missing arguments in `FcnPrepare!` are added;
  - the column kernels in `Solve!` no longer exceed 1024 threads for Nz > 102;
  - `CGPARAM="name=value,..."` overrides problem parameters.

Column tests (`Examples/PBLRun.jl`, `Jobs/NHCart/CG/ColumnDiffusion`, `ColumnTKE`):

| Test | Result |
|---|---|
| Gaussian u bump, K = 10 m²/s, 1 h, Δz = 20 m | implicit Δt = 60 s: peak 1.7462 against 1.7461 analytic; same error as explicit Δt = 5 s |
| θ diffusion, N = 0.02 1/s | implicit 60 s against explicit 5 s: 5·10⁻⁵ K |
| neutral, 10 m/s, z0 = 0.1 m, TKE, 2 h | u\* = κU₁/ln(z₁/z0); log law within 1–2 % up to 90 m; K = 0.85–0.9 κu\*z; e = 0.85–1.0 u\*²/c_m²; implicit 60 s against explicit 1 s: 10⁻³ |
| heated / cooled surface (±5 K) | unstable: e up to 0.9 m²/s², K about 30 m²/s; stable: shallow layer, K about 1 m²/s |
| Float32 against Float64 (unstable column) | 10⁻⁴ |

Next: the Doyle–Durran rotor case (`Jobs/NHCart/CG/RotorCart`, periodic in x), then the
DGSEM steps of section 4.

## 1 What exists

| Part | Where | Used by |
|---|---|---|
| Switches `--SurfaceFlux`, `--SurfaceFluxMom`, `--VerticalDiffusion`, `--VerticalDiffusionMom`, `--JacVerticalDiffusion`, `--Turbulence`, `--TkePos`, `--SurfaceScheme` | `src/Parameters/parse_commandline.jl`, set up in `Examples/DriverCG.jl` (about lines 505–550) | CGSEM; `DriverDG.jl` reads them but DGSEM never uses them |
| Surface values (surface temperature, saturation vapour density) | `src/Surfaces/surface.jl`: `DefaultSurface`, `HeldSuarezMoistSurface` | moist Held–Suarez, moist baroclinic wave, Frierson jobs |
| Surface exchange coefficients | `MOSurfaceFlux` with Businger functions (`functions.jl`) and a Newton iteration for the stability parameter (`solver.jl`, `MOSTIteration`, `MOSTSeaIteration`); `HeldSuarezMoistSurfaceFlux` with constant coefficients | `--SurfaceScheme="MOST"` |
| Land use | `LandUseData`: two classes; z0M = z0H = 0.01 m; every column is set to sea (`LandClass = 2`) | all |
| Surface flux tendencies | `Surfaces.SurfaceFlux` (ρθ, and ρ and ρqv for moist runs), `SurfaceFluxScalarsKernel!`, `SurfaceFluxMomentumKernel!` (`CGSEM/MomentumKernel.jl`) | explicit |
| Eddy diffusivity | `SimpleKoefficient` (Reed–Jablonowski-type K = C_E·u·Δz/2, decaying above p_pbl), `TkeKoefficient` (`Surfaces/BoundaryLayer.jl`), `EddyCoefficientKernel!` (`CGSEM/DiagnosticKernel.jl`) | |
| Vertical diffusion | `VerticalDiffusionScalarKernel!` (ρθ, tracers, TKE; `OperatorKernel.jl`), `VerticalDiffusionMomentumKernel!` (u, v) | scalars implicit with `--JacVerticalDiffusion`, momentum explicit |
| Implicit diffusion | `JacTransportKernel!` builds a tridiagonal `JDiff`; `SchurSolveTriDiffKernel!` solves it for the listed scalars before the acoustic Schur solve (`CGSEM/Solve.jl`) | Rosenbrock W-method |
| TKE | `Models.TKEModel` (`src/Models/Turbulence.jl`): shear production, Richardson-number correction, dissipation; `TkeSourceKernel!` adds the source and writes K at the interfaces | `Jobs/NHCart/CG/BoxTKE`, `TkeCart` |
| Not wired in or not working | `FriersonSurfaceFlux`, `FriersonKoefficient`, `BLHRichardson` (undefined variables); `SeaSurfaceDrag.jl` and `wrfscheme.jl` (half-translated Fortran; `SeaSurfaceDrag.jl` is not included) | – |

Jobs using the physics: `NHSphere/CG/HeldSuarezMoistSphere*` and
`BaroWaveMoistSphere_32Elem` (MOST, scalar diffusion with the simple K, no
TKE), `FriersonSphere_32Elem`, and `NHCart/CG/BoxTKE` and `TkeCart` (MOST and
TKE). According to `Jobs/README.md` they are only checked for three time
steps with a finite solution. None has been validated against a reference.

## 2 Problems found

These come from reading the code. Each must be confirmed by a test before it is fixed.

1. **Bulk formula with the wrong velocity.** For MOST, `MOSurfaceFlux` stores
   the friction velocity u\* in `SD[uStarPos]`. `SurfaceFlux` then computes
   the heat and moisture flux as `CT * uStar * (T - TS)`. But `CT` is a bulk
   coefficient that belongs to the wind speed |v|. With u\* the flux is too
   small by about u\*/|v| = √C_M, roughly a factor of 20. `SimpleKoefficient`
   also uses u\* where Reed–Jablonowski use |v|. With the constant-coefficient
   scheme (`HeldSuarezMoistSurfaceFlux`) the slot holds |v|, so there the
   formula is right. The momentum kernel computes |v| itself, so it is consistent.
2. **Reference height of MOST.** `SurfaceFluxData!` passes `dz[1,IC]`, the
   thickness of the lowest layer, as z. The state of the lowest layer sits at
   its middle, z = Δz₁/2.
3. **Factor 2 in the momentum surface flux.** `SurfaceFluxMomentumKernel!`
   divides by Δz₁/2. The scalar flux divides by Δz₁. One of the two is
   inconsistent.
4. **Wrong K for momentum diffusion.** `Cache.KV` has the shape (M, Nz, NumG).
   `VerticalDiffusionMomentumKernel!` reads it as `K[Iz,IC]`, i.e. with linear
   indexing, while the scalar kernel reads `K[1,Iz,IC]` correctly.
5. **K at layer midpoints used as an interface value.** `EddyCoefficientKernel!`
   computes K per layer (Iz). The diffusion kernels use it as the value at the
   interface between Iz and Iz+1, which is half a layer off.
6. **Diffusion Jacobian.** In `JacTransportKernel!` the lower density `RhoB`
   is read from `Rho[1,Iz+1,IC]` (should be `Iz`). The off-diagonal entries
   should be checked against the explicit operator.
7. **Surface data called with the wrong arguments.** In `FcnPrepare!`, the branch
   "vertical diffusion without surface fluxes" calls `SurfaceFluxData!` with
   arguments missing (no `T`, `PotT`). It would fail at run time.
8. **TKE model.**
   - Buoyancy is switched off (`N2 = 0.0` after computing it).
   - The stable length scale uses √e/N² instead of √(e/N²).
   - The dissipation coefficient `CDis` mixes two formulations; `Diss` is
     computed but not used.
   - `LenScale = min(Δz, 200 m)` makes the mixing length depend on the grid.
   - There is no surface boundary condition for e (u\*-based) and no
     dissipation-limited time-step treatment.
   - In `TkeSourceKernel!` the top cell gets only half of the interface
     source, because the branch `Iz == Nz` can never run.
   - TKE is not kept positive.
9. **Float64 literals in Float32 runs.** `max(uT, 0.1)` (`UzCoefficientGPU`),
   `0.5 *` in `SurfaceFluxMomentumKernel!`, `LenScale = 100.0`, `N2 = 0.0` and
   `max(... ,1.e-1)` promote to Float64. The results stay correct, but on the
   GPU they cost speed and registers (compare the kernel work in BaroWave10Days).
10. **Stability of explicit momentum diffusion.** With 10–20 m layers near the
    ground and K of 10–100 m²/s, explicit diffusion is limited to
    Δt < Δz²/(2K), i.e. a few seconds. The dynamics runs at 150 s. Momentum
    diffusion and surface drag must therefore be implicit, like the scalars.

## 3 Plan for CGSEM

**Step 1: single-column test bed.** Build the case on the existing
`VerticalColumn` setup, with a doubly periodic small box or one column.
Compare against analytic solutions or published single-column results:
- neutral Ekman layer with constant K (analytic solution);
- GABLS1 (stable, prescribed cooling of the surface);
- a dry convective boundary layer with a prescribed surface heat flux.

The same cases later serve DGSEM. Outputs: profiles of u, v, θ and K,
the surface fluxes, and the boundary-layer height.

**Step 2: fixes (items 1–9).**
- One place computes the surface exchange: |v| and u\* separately,
  z = Δz₁/2, and C_M, C_H from MOST.
- The fluxes enter the lowest cell as the lower boundary condition of the
  diffusion operator, rather than as a separate source with its own factor.
- K lives at the interfaces, like w, and the kernels index it consistently.
- Float32-clean literals.

Check: the Ekman layer with constant K must converge to the analytic
solution, and the surface flux must match a hand calculation.

**Step 3: implicit momentum diffusion and drag.**
- Add tridiagonal solves for u and v with the same structure as `JDiff`
  (`ρ`-weighted, per component).
- Put the linearised surface drag C_M|v|/Δz₁ on the diagonal of the lowest row.
- Like the scalars, they are solved before the acoustic Schur solve, an
  approximate factorisation of the W-matrix.

Check: stable at Δt = 150 s with K up to 100 m²/s and Δz₁ = 20 m, and the
same result as an explicit run with small Δt.

**Step 4: TKE closure.** A 1.5-order closure with a prognostic e at the
interfaces (staggered like w, so shear and N² are available there without
averaging):
- production K_m S² − K_h N²;
- dissipation c_ε e^{3/2}/l;
- K_m = c_m l √e and K_h = K_m / Pr_t(Ri);
- mixing length: Blackadar, limited by 0.76√(e/N²) in stable layers, or
  Bougeault–Lacarrère;
- surface value e = u\*²/√(c_m c_ε);
- dissipation implicit (linearised), production explicit;
- TKE vertical diffusion through the same tridiagonal solve;
- clipping or a positivity-preserving update of e.

Which closure exactly — ICON TURBDIFF-like, Mellor–Yamada level 2.5 or
Deardorff/Bougeault–Lacarrère — is a decision for you (see section 6).

Check: GABLS1 and the convective boundary layer against published results;
then `TkeCart` (flow over a hill).

**Step 5: 3D tests on the sphere.**
- Held–Suarez with real surface fluxes and diffusion instead of the Rayleigh
  friction: the Reed–Jablonowski simple physics for the moist Held–Suarez and
  moist baroclinic wave cases, where results exist for comparison.
- Then the moist baroclinic wave with TKE.
- Measure the cost per step, in Float32 and Float64.

## 4 Bringing it to DGSEM

**Principle.** The pointwise physics — MOST, the bulk fluxes, the TKE sources,
K from e and l — stays in shared functions, as the flux functions are now
selected through `FluxSelect`. Only the discretisation of the vertical
operator and the boundary condition is DGSEM-specific. The column test bed
of step 1 runs both dycores.

**A. Vertical resolution first.** The DGSEM jobs use vertical elements of
2.31 km with 5 GLL nodes. The first node above the ground is about 400 m
high. A boundary layer needs a stretched vertical element grid, with a lowest
element of about 100–200 m and therefore a first node gap of 20–35 m.

Check that DGSEM runs with stretched vertical elements, and what that costs:
- the time step is horizontal, so it is unchanged;
- the vertical acoustics are implicit (HDG), so they are not affected.

**B. Surface fluxes as a boundary flux.** DGSEM has a node on the ground, and
the bottom face already gets a numerical flux (the vertical Riemann solver
with a wall).
- The surface momentum, heat and moisture fluxes are added to this
  boundary flux and lifted to the bottom node. That is the natural DG weak
  form; no source term with a Δz factor is needed.
- For the momentum flux, the Cartesian velocity is projected onto the local
  tangent plane with the surface normal `nS`, the same projection CGSEM uses.
- MOST needs values at a height above the ground. Either use the first node
  above the ground, or evaluate the element polynomial at a fixed reference
  height (e.g. 10 m). The second choice is independent of the vertical grid
  and is recommended.

**C. Vertical diffusion with HDG.** Diffusion is the textbook case for HDG.
For a scalar φ, write the problem as q = K ∂zφ, ∂t(ρφ) = ∂z(ρq), with a
single trace value of φ at each element interface and a stabilisation
parameter τ ∼ K/Δz.
- Static condensation gives a tridiagonal system per column for the traces,
  the same structure the acoustic HDG solve already has.
- First version: a separate implicit diffusion solve per variable (θ,
  tracers, the two tangential velocity components, e) before the acoustic
  solve, like the CGSEM approximate factorisation. It reuses the fused
  column-kernel pattern of the HDG solve.
- Later: coupling diffusion and acoustics in one W-matrix, if splitting
  errors show up.
- Explicit alternative for testing only: BR1/LDG with an auxiliary gradient,
  limited by the diffusion time step.

**D. TKE in DGSEM.**
- Start with a diagnostic closure: Mellor–Yamada level 2, where K follows from
  the local shear, N² and the mixing length, with no prognostic e. It needs
  no new variable and no positivity treatment, and it tests C and B.
- Then a prognostic e as an additional tracer at the nodes. DG solutions
  oscillate, so e needs a positivity-preserving limiter (Zhang–Shu type) or
  a transformed variable. The shear and N² come from the element
  derivative matrix, which DGSEM already has.

**E. Equation sets.** In the vector-invariant DGSEM the potential temperature
is diffused through ρθ, consistently with θ = ρθ/ρ. The conservative forms
(Cons, Diff) are handled the same way. Diffusing the Cartesian momentum must
not act on the vertical velocity; only the tangential part is diffused.

**F. Tests.**
- The same column cases as for CGSEM, using `Jobs/NHCart/DG/VerticalColumnDG`.
- Then `HeldSuarezDrySphereDGQuad_20Elem` with surface fluxes and diffusion
  instead of the Rayleigh friction, compared with CGSEM on the same case.

## 5 Order of work

| # | Work | Dycore | Depends on |
|---|---|---|---|
| 1 | Column test bed (Ekman, GABLS1, convective BL) | CGSEM | – |
| 2 | Fixes of section 2, Float32 clean-up | CGSEM | 1 |
| 3 | Implicit momentum diffusion and drag | CGSEM | 2 |
| 4 | TKE closure | CGSEM | 3 |
| 5 | Shared pointwise physics module | both | 2 |
| 6 | Stretched vertical DG elements | DGSEM | – (in parallel with 1–4) |
| 7 | Surface fluxes as a boundary flux | DGSEM | 5, 6 |
| 8 | HDG vertical diffusion | DGSEM | 7 |
| 9 | Diagnostic, then prognostic TKE | DGSEM | 4, 8 |
| 10 | Held–Suarez / moist baroclinic wave with physics, CGSEM against DGSEM, timings | both | 4, 9 |

## 6 Decisions for you

- **TKE closure:** ICON TURBDIFF-like (Raschendorfer), Mellor–Yamada 2.5, or a
  simpler Deardorff/Bougeault–Lacarrère scheme. This decides how close the
  results can be compared with ICON.
- **Reference cases:** which cases and data to compare against — GABLS, DYCOMS,
  Reed–Jablonowski simple physics, ICON single-column runs.
- **Land surface:** whether more land-use classes and a land surface are
  needed, or prescribed surface temperature and roughness are enough for
  now.
- **EDMF:** CGSEM has EDMF variables (`wEDMF`, `aRhoEDMF`). They were not
  reviewed here. Should the TKE closure later be the environment part of an
  EDMF scheme?
