function AddVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
Grid.dzeta = zeros(nz)
Grid.H = H
@. Grid.dzeta = H/nz;
Grid.H=H
for i=2:nz+1
  Grid.z[i] = Grid.z[i-1] + Grid.dzeta[i-1]
end
for i=1:nz
  Grid.dzeta[i] = Grid.z[i+1] - Grid.z[i]  
  Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
end
end

function AddStretchICONVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64,sigma::Float64,lambda::Float64)
  Grid.zP=zeros(nz);
  Grid.z=zeros(nz+1);
  Grid.dzeta = zeros(nz)
  Grid.H = H
  for iz=1:nz
    i = nz + 1 - iz  
    Grid.z[iz+1] = H*(2.0/pi*acos(((i - 1) / nz)^sigma))^lambda
  end
# Grid.z[1] = 0.0
# Grid.z[2] = 20.0
# dz = (H-Grid.z[2]) / (nz - 1)
# for i = 2 : nz -1 
#   Grid.z[i+1] = Grid.z[i] + dz  
# end
  Grid.z[nz+1] = H
  for i=1:nz
    Grid.dzeta[i] = Grid.z[i+1] -Grid.z[i]
    Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
  end
end

function AddStretchDCMIPVerticalGrid!(Grid::GridStruct, nz::Int,H::Float64, mue::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
Grid.dzeta = zeros(nz)
Grid.H = H
for iz=1:nz
  i = nz + 1 - iz  
  Grid.z[iz+1] = H * (sqrt(mue * ((nz - i +1) / nz)^2 + 1.0) - 1.0) / (sqrt(mue + 1.0) - 1.0)
end
for i=1:nz
  Grid.dzeta[i] = Grid.z[i+1] -Grid.z[i]
  Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
end
end

function AddExpStretchVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64,
  dz_bottom::Float64,dz_top::Float64)

  Grid.zP=zeros(nz);
  Grid.z=zeros(nz+1);
  Grid.dzeta = zeros(nz)
  Grid.H = H

  z_bottom = 0.0
  z_top = H
  nelems = nz

  # define the inverse σ⁻¹ exponential stretching function
  exp_stretch(ζ, h) = ζ == 1 ? ζ : -h * log(1 - (1 - exp(-1 / h)) * ζ)

  # nondimensional vertical coordinate (]0.0, 1.0])
  ζ_n = LinRange(1.0, nelems, nelems) / nelems

  # find bottom height variation
  find_bottom(h) = dz_bottom - z_top * exp_stretch(ζ_n[1], h)
  # we use linearization
  # h_bottom ≈ -dz_bottom / (z_top - z_bottom) / log(1 - 1/nelems)
  # to approx bracket the lower / upper bounds of root sol
  guess₋ = -dz_bottom / (z_top - z_bottom) / log(1 - Float64(1 / (nelems - 1)))
  guess₊ = -dz_bottom / (z_top - z_bottom) / log(1 - Float64(1 / (nelems + 1)))
  h_bottom_sol = RootSolvers.find_zero(
        find_bottom,
        RootSolvers.SecantMethod(guess₋, guess₊),
        RootSolvers.CompactSolution(),
        RootSolvers.ResidualTolerance(1e-3),
    )
  if h_bottom_sol.converged !== true
    error(
        "h_bottom root failed to converge for dz_bottom: $dz_bottom on domain ($z_bottom, $z_top)",
    )
  end
  h_bottom = h_bottom_sol.root

  # find top height variation
  find_top(h) = dz_top - z_top * (1 - exp_stretch(ζ_n[end - 1], h))
  # we use the linearization
  # h_top ≈ (z_top - dz_top) / z_top / log(nelem)
  # to approx braket the lower, upper bounds of root sol
  guess₋ = ((z_top - z_bottom) - dz_top) / z_top / Float64(log(nelems + 1))
  guess₊ = ((z_top - z_bottom) - dz_top) / z_top / Float64(log(nelems - 1))
  h_top_sol = RootSolvers.find_zero(
      find_top,
      RootSolvers.SecantMethod(guess₋, guess₊),
      RootSolvers.CompactSolution(),
      RootSolvers.ResidualTolerance(Float64(1e-3)),
  )
  if h_top_sol.converged !== true
      error(
          "h_top root failed to converge for dz_top: $dz_top on domain ($z_bottom, $z_top)",
      )
  end
  h_top = h_top_sol.root

  # scale height variation with height
  h =
    h_bottom .+
    (ζ_n .- ζ_n[1]) * (h_top - h_bottom) / (ζ_n[end - 1] - ζ_n[1])
    faces = (z_bottom + (z_top - z_bottom)) * exp_stretch.(ζ_n, h)

  # add the bottom level
  faces = [z_bottom; faces...]
  for iz = 1 : nz + 1
    Grid.z[iz] = faces[iz]
  end
  for i = 1 : nz
    Grid.dzeta[i] = Grid.z[i+1] - Grid.z[i]
    Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
  end
end
