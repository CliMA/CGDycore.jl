# Jobs

Each file is a bash script that starts one example with all its options.
Run it from the repository root:

```bash
export JuliaDevice="GPU" JuliaGPU="CUDA"   # or JuliaGPU="Metal"; default CPU
bash Jobs/NHSphere/DG/BaroWaveDrySphereDGQuad_VI
```

Options on the command line after the job's options override them
(ArgParse takes the last value), e.g. `--FloatTypeBackend="Float32"`.

## Layout

| Directory | Model | Dycore | Driver |
|---|---|---|---|
| `NHCart/CG` | nonhydrostatic, Cartesian box | CGSEM | `Examples/DriverCG.jl` |
| `NHCart/DG` | nonhydrostatic, Cartesian box | DGSEM | `Examples/DriverDG.jl` |
| `NHSphere/CG` | nonhydrostatic, sphere | CGSEM | `Examples/DriverCG.jl` |
| `NHSphere/DG` | nonhydrostatic, sphere | DGSEM | `Examples/DriverDG.jl` |
| `SWSphere/CG` | shallow water, sphere | CGSEM | `Examples/DriverCG.jl` |
| `FEM` | shallow water, mixed finite elements | FEM | `Examples/FEM/...` |
| `FV` | finite volumes | FV | `Examples/FV/...` |
| `Tools` | Earth topography on the grid | | `Examples/CreateOrography.jl` |

Cartesian CGSEM jobs need `--GridForm="Cartesian"`, otherwise the driver
builds the cubed sphere.

## Notes

- Earth orography (`--TopoS="EarthOrography"`) reads
  `ETOPO_2022_v1_60s_N90W180_surface.nc` from the working directory.
- `NHSphere/CG/HeldSuarezMoistSphereAMD` starts through the LUMI
  `gpu_wrapper` and runs only there.
- DGSEM needs the pressure and the geopotential in the auxiliary fields:
  `--pAuxPos=1 --GPAuxPos=2`.
- `--Forcing=true` (Held-Suarez) needs `--pAuxPos=1`, surface fluxes or
  vertical diffusion need `--SurfaceScheme="MOST"`.

## Status

All CGSEM and DGSEM jobs run three time steps on an NVIDIA A100 (Levante)
with a finite solution; the two Earth orography jobs
(`NHSphere/CG/BaroWaveDrySphereOro`, `NHSphere/CG/HeldSuarezMoistSphereOro`)
run one simulated hour. Not tested: `NHSphere/CG/HeldSuarezMoistSphereAMD`
(LUMI only) and the FEM and FV jobs.
