function AddVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
Grid.H = H
@. Grid.dzeta = H/nz;
Grid.H=H
for i=2:nz+1
  Grid.z[i] = Grid.z[i-1] + Grid.dzeta[i-1]
end
for i=1:nz
  Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
end
end

function AddStretchICONVerticalGrid!(Grid::GridStruct,nz::Int,H::Float64,sigma::Float64,lambda::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
Grid.H = H
for iz=1:nz
  i = nz + 1 - iz  
  Grid.z[iz+1] = H*(2.0/pi*acos(((i - 1) / nz)^sigma))^lambda
end
for i=1:nz
  Grid.dzeta[i] = Grid.z[i+1] -Grid.z[i]
  Grid.zP[i] = 0.5 * (Grid.z[i] + Grid.z[i+1])
end
end

function AddStretchDCMIPVerticalGrid!(Grid::GridStruct, nz::Int,H::Float64, mue::Float64)
Grid.zP=zeros(nz);
Grid.z=zeros(nz+1);
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

#=
function IntervalMesh(
    domain::IntervalDomain{CT},
    stretch::ExponentialStretching{FT};
    nelems,
) where {CT <: Geometry.Abstract1DPoint{FT}} where {FT}
    if nelems < 1
        throw(ArgumentError("`nelems` must be ≥ 1"))
    end
    cmin = Geometry.component(domain.coord_min, 1)
    cmax = Geometry.component(domain.coord_max, 1)
    R = cmax - cmin
    h = stretch.H / R
    η(ζ) = -h * log1p((expm1(-1 / h)) * ζ)
    faces =
        [CT(cmin + R * η(ζ)) for ζ in range(FT(0), FT(1); length = nelems + 1)]
    IntervalMesh(domain, faces)
end
=#
