@inline function VelSphere2Cart(VelSp,lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = eltype(lon)(0)
  rot23 = coslat
  rot33 = sinlat
# VelCa = rot'*VelSp
  VelCa1 = rot11 * VelSp[1] + rot21 * VelSp[2] + rot31 * VelSp[3]
  VelCa2 = rot12 * VelSp[1] + rot22 * VelSp[2] + rot32 * VelSp[3]
  VelCa3 = rot13 * VelSp[1] + rot23 * VelSp[2] + rot33 * VelSp[3]
  return SVector{3}(VelCa1, VelCa2, VelCa3)
end

@inline function VelCart2Sphere(VelCa,lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = eltype(lon)(0)
  rot23 = coslat
  rot33 = sinlat
# VelSp = rot*VelCa
  VelSp1 = rot11 * VelCa[1] + rot12 * VelCa[2] + rot13 * VelCa[3]
  VelSp2 = rot21 * VelCa[1] + rot22 * VelCa[2] + rot23 * VelCa[3]
  VelSp3 = rot31 * VelCa[1] + rot32 * VelCa[2] + rot33 * VelCa[3]
  return SVector{3}(VelSp1, VelSp2, VelSp3)
end

@inline function MSphere2Cart(lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = eltype(lon)(0)
  rot23 = coslat
  rot33 = sinlat
  return  @SMatrix [rot11 rot21 rot31; rot12 rot22 rot32; rot13 rot23 rot33]
end

@inline function MCart2Sphere(lon,lat)
  sinlon = sin(lon)
  coslon = cos(lon)
  sinlat = sin(lat)
  coslat = cos(lat)
  rot11 = -sinlon
  rot21 = -sinlat * coslon
  rot31 = coslat * coslon
  rot12 = coslon
  rot22 = -sinlat * sinlon
  rot32 = coslat * sinlon
  rot13 = eltype(lon)(0)
  rot23 = coslat
  rot33 = sinlat
  return  @SMatrix [rot11 rot12 rot13; rot21 rot22 rot23; rot31 rot32 rot33]
end

@kernel function spherical_to_cartesian_velocity_kernel!(
    v_x, v_y, v_z,                       # Outputs (M, nz, NDoF)
    v_r, v_θ, v_φ,                       # Inputs  (M, nz, NDoF)
    θ_arr, φ_arr                         # Inputs (NDoF,)
)
    idx = @index(Global)
    
    if idx <= size(v_x, 3)
        @fastmath begin
            # Berechne sin und cos simultan direkt in die Register (Null Speicher-Traffic)
            s_θ, c_θ = sincos(θ_arr[idx])
            s_φ, c_φ = sincos(φ_arr[idx])
            
            # Da M (<=6) und nz (<=40) klein sind, loopen wir linear intern.
            for nz_idx in 1:size(v_x, 2)
                for m_idx in 1:size(v_x, 1)
                    
                    # Nur noch die 3 Geschwindigkeitskomponenten laden
                    dot_r = v_r[m_idx, nz_idx, idx]
                    dot_θ = v_θ[m_idx, nz_idx, idx]
                    dot_φ = v_φ[m_idx, nz_idx, idx]
                    
                    # Vereinfachte Terme für r = 1
                    dot_θ_sin_φ = s_φ * dot_θ
                    dot_φ_cos_φ = c_φ * dot_φ
                    dot_r_sin_φ = dot_r * s_φ
                    
                    # Transformation (wird vom Compiler perfekt zu FMA-Instruktionen reduziert)
                    v_x[m_idx, nz_idx, idx] = dot_r_sin_φ * c_θ - dot_θ_sin_φ * s_θ + dot_φ_cos_φ * c_θ
                    v_y[m_idx, nz_idx, idx] = dot_r_sin_φ * s_θ + dot_θ_sin_φ * c_θ + dot_φ_cos_φ * s_θ
                    v_z[m_idx, nz_idx, idx] = dot_r * c_φ - s_φ * dot_φ
                end
            end
        end
    end
end

using KernelAbstractions

@kernel function cartesian_to_spherical_velocity_transposed_kernel!(
    v_r, v_θ, v_φ,                       # Outputs (M, nz, NDoF) - lineare Komponenten
    v_x, v_y, v_z,                       # Inputs  (M, nz, NDoF)
    θ_arr, φ_arr                         # Inputs  (NDoF,)
)
    idx = @index(Global)

    if idx <= size(v_r, 3)
        @fastmath begin
            # Berechne sin und cos on-the-fly direkt in die CPU/GPU-Register
            s_θ, c_θ = sincos(θ_arr[idx])
            s_φ, c_φ = sincos(φ_arr[idx])

            # Da M (<=6) und nz (<=40) klein sind, loopen wir linear intern.
            for nz_idx in 1:size(v_r, 2)
                for m_idx in 1:size(v_r, 1)

                    # Kartesische Geschwindigkeiten laden
                    dot_x = v_x[m_idx, nz_idx, idx]
                    dot_y = v_y[m_idx, nz_idx, idx]
                    dot_z = v_z[m_idx, nz_idx, idx]

                    # Inverse Transformation via exakt transponierter Matrix
                    # Zeile 1: dot_r = sin(φ)cos(θ)*ẋ + sin(φ)sin(θ)*ẏ + cos(φ)*ż
                    v_r[m_idx, nz_idx, idx] = (s_φ * c_θ) * dot_x + (s_φ * s_θ) * dot_y + c_φ * dot_z

                    # Zeile 2: dot_θ = -sin(θ)*ẋ + cos(θ)*ẏ
                    v_θ[m_idx, nz_idx, idx] = -s_θ * dot_x + c_θ * dot_y

                    # Zeile 3: dot_φ = cos(φ)cos(θ)*ẋ + cos(φ)sin(θ)*ẏ - sin(φ)*ż
                    v_φ[m_idx, nz_idx, idx] = (c_φ * c_θ) * dot_x + (c_φ * s_θ) * dot_y - s_φ * dot_z
                end
            end
        end
    end
end
:
