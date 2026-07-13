using MuladdMacro
using Trixi: get_node_vars, set_node_vars!, get_node_coords, @threaded

"""
    JacobianCacheAndres

    Data structure to store the Jacobian matrices for DGSEM-HEVI discretizations of the  
    CompressibleEulerPotentialTemperatureEquationsWithGravity2D and
    CompressibleEulerPotentialTemperatureEquationsWithGravity3D.
    
    This Jacobian structure applies the following permutation to a discretization that uses
    `nz` elements of polynomial degree `polydeg` in the vertical direction:

    The solution array is permutated, such that it contains:
    * The density DOFs of 1 to polydeg of element 1, then
    * The potential temperature DOFs 1 to polydeg of element 1, then 
    * The vertical momentum DOFs 1 to polydeg of element 1, then 
    * ...
    * The density DOFs of 1 to polydeg of element `nz`, then
    * The potential temperature DOFs 1 to polydeg of element `nz`, then 
    * The vertical momentum DOFs 1 to polydeg of element `nz`, then 
    * The density DOFs at node polydeg+1 of element 1
    * The vertical momentum DOFs DOFs at node polydeg +1 element 1
    * The potential temperature DOFs at node polydeg+1 element 1
    * ...
    * The density DOFs at node polydeg+1 element `nz`
    * The vertical momentum DOFs DOFs at node polydeg +1 element `nz`
    * The potential temperature DOFs at node polydeg+1 element `nz`

    This leads to the following sparsity pattern  (example for `nz` = 4, polydeg = 5):

    ⎡⠑⢄⠀⠀⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⢸⠀⠀⠀⠀⠀⎤
    ⎢⠀⠀⠑⢄⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⢸⠀⠀⠀⠀⠀⎥
    ⎢⢄⠀⢠⣤⣵⢟⠛⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⠘⡄⠀⠀⠀⠀⎥
    ⎢⠀⠑⠼⠿⠿⠀⠑⢄⠀⠀⡀⠀⢀⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⢀⡇⡀⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⢸⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⠀⠀⡇⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⢸⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⠈⠁⡇⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠐⢄⠀⣶⣶⡟⢍⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⠐⠂⢱⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠛⠛⠃⠀⠑⢄⠀⠠⠀⠀⣤⣤⡄⠀⠀⠀⠀⠀⠀⠀⎥⠀⠀⠼⢠⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⎥⠀⠀⣀⢸⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⎥⠀⠀⠀⢸⠀⠀⎥  = ⎡A B⎤
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⢸⣿⣿⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⎥⠀⠀⠉⠀⡇⠀⎥    ⎣C D⎦
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠉⠀⠀⠑⢄⠀⠂⠀⢰⣶⣶⎥⠀⠀⠀⠐⠃⡆⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⢸⣿⣿⎥⠀⠀⠀⠠⠄⡇⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⣀⣑⣼⠿⠿⎥⠀⠀⠀⢀⡀⢇⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⣿⣿⡇⠑⢄⎥⠀⠀⠀⠀⠀⢸⎥
    ⎢------------------------------⎥------⎥
    ⎢⠀⠀⠐⠒⠒⠭⠭⠅⠀⠀⠇⠀⢸⣀⣀⠀⠀⢀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥⠻⢇⣀⠀⠀⠀⎥
    ⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠑⠒⠒⠀⠀⢘⣀⣀⠧⠤⠄⠀⠀⡄⠀⢠⠀⠀⎥⠀⠈⠛⣤⡄⠀⎥
    ⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠉⠁⠀⠀⠥⠤⢜⣒⣒⎥⠀⠀⠀⠈⠱⣦⎦

    ### A-block structure

    A is itself a 3x3 block matrix and we write as:
    ⠀⢄⠀⠠⠀⠀⣤⣤⡄
     ⠀⠑⢄⠀⠀⣿⣿⡇   A11 A12 A13
     ⠀⠀ ⠑⢄⣿⣿⡇ = A21 A22 A33
     ⠑⢄⢸⣿⣿⠑⢄⠀   A31 A32 A33 
     ⠀⠀⠉⠉⠉⠀⠀⠂
    
    A11 = fac·I             (ρ–ρ,   diagonal, uniform)
    A12 = A12_11·e₁e₁ᵀ     (ρ–ρθ,  sparse: only [1,1] nonzero)
    A13                     (ρ–ρw,  full matrix)
    A21 = 0
    A22 = diag(A22_diag)   (ρθ–ρθ, diagonal, non-uniform)
    A23                     (ρθ–ρw, full matrix)
    A31 = FacGrav·I        (ρw–ρ,  diagonal = gravity in momentum eq)
    A32                     (ρw–ρθ, full matrix)
    A33 = diag(A33_diag)   (ρw–ρw, diagonal, non-uniform; absorbed into SA)

    We keep the same 3x3 block structure for B, C, D.

    ### Row orderings
    A matrix rows: (ρ [polydeg entries], ρθ [polydeg entries], ρw [polydeg entries]) for each element   →  variable 1 = ρ, 2 = ρθ, 3 = ρw
    C matrix rows: (ρ, ρw, ρθ)   →  edge eq  1 = ρ, 2 = ρw, 3 = ρθ

    The two orderings differ only in positions 2 and 3 (ρθ ↔ ρw).
    Changing C to match A would be pure row/column relabelling in B, C, D
    with identical nonzero count and FLOPs — no performance gain.

    ### B-block naming  (row group = A-variable index)
    B1 = ρ  interior block (A-row group 1)
    B2 = ρθ interior block (A-row group 2)
    B3 = ρw interior block (A-row group 3)
    B1_col2, B2_col2: ρ and ρθ interior driven by ρw-edge DOF (col 2)
    B3_col3:          ρw interior driven by ρθ-edge DOF (col 3)

    ### C-block naming  (col group = A-variable index; row = C-row-order index)
    C3_rows13[1,j] : C[ρ-edge  eq (row 1), ρw-interior j]  (col group 3 = ρw)
    C3_rows13[2,j] : C[ρθ-edge eq (row 3), ρw-interior j]
    C2_row2[j]     : C[ρw-edge eq (row 2), ρθ-interior j]  (col group 2 = ρθ)
    FacGrav·r1_1   : C[ρw-edge eq (row 2), ρ-interior node 1]  (gravity in C)
    
"""
mutable struct JacobianCacheAndres{NZ, M, RealT, AType2 <: AbstractArray, AType3 <: AbstractArray,
                             AType4 <: AbstractArray, zColType <: AbstractArray,
                             JacobianAuxType,
                             RotationMatrixType, IDMapType, IZMapType, UVertType}
    M::Int # TODO This is M in Oswald's code and can be substituted with eachnode
    nz::Int
    nx::Int
    fac::RealT
    A12_11::AType2
    A22_diag::AType3
    A33_diag::AType3
    FacGrav::RealT
    A13::AType4
    A23::AType4
    A32::AType4
    B1m_row1_23::AType3
    B2m_row1_23::AType3
    B3m_row1_23::AType3
    B1_col2::AType3
    B2_col2::AType3
    B3_col3::AType3
    C2_row2::AType3
    C3_rows13::AType4
    C2p_col1::AType3
    C3p_col1::AType3
    SA::AType4
    # Storage for the Schur complement AND the D matrix of the vertical Jacobians.
    # We use a block-tridiagonal matrix format with the known sparsity pattern: 
    #   * Diagonal blocks have D[i][3,1] = 0
    #   * Off-diagonal blocks L[i], U[i]:  first column is always zero.)
    SchurD  :: AType4   # (3,3,nz,NumG)     # Diagonal bloks (3x3)
    SchurL  :: AType4   # (3,2,nz-1,NumG)   # Upper blocks (3x2)
    SchurU  :: AType4   # (3,2,nz-1,NumG)   # Lower blocks (3x2)
    # Fields for LU facotization and solve... See block_lu!(F::JacobianCacheAndres)
    Schurm      :: AType2      # (nz,NumG)
    Schurd1inv  :: AType2      # (nz,NumG)
    Schurd1vw   :: AType3     # (2,nz,NumG)
    Schuru1     :: AType3     # (2,nz-1,NumG)
    Schurl1     :: AType3     # (2,nz-1,NumG)
    SchurLhat   :: AType4     # (2,2,nz-1,NumG)
    SchurDinv   :: AType4     # (2,2,nz,NumG)
    SchurUhat   :: AType4     # (2,2,nz-1,NumG)
    #
    rs      :: AType2
    NumG::Int
    dz::zColType
    cS::RealT
    jacobian_aux::JacobianAuxType
    wPos::Int
    ThPos::Int
    rotation_matrix::RotationMatrixType
    ID_map::IDMapType
    iz_map::IZMapType
    u_vert::UVertType
end

# Constructor
function JacobianCacheAndres(M, nz, nx, NumG, semi, cS, dt, 
                             equations::Union{CompressibleEulerPotentialTemperatureEquationsWithGravity2D, 
                             CompressibleEulerPotentialTemperatureEquationsWithGravity3D})
    FT = Float64
    polydeg = M - 1
    M2 = M - 2 # To be deprecated
    fac = inv(dt) # Value of MOST entries in diagonal of A11, A22, and A33 (diagonal matrices)
    A12_11 = zeros(FT, nz, NumG) # A12[1,1]
    A22_diag = zeros(FT, polydeg, nz, NumG) # diagonal of A22
    A33_diag = zeros(FT, polydeg, nz, NumG) # diagonal of A33
    A22_diag .= fac # All entries are fac but the first, which we will overwrite
    A33_diag .= fac # All entries are fac but the first, which we will overwrite
    FacGrav = semi.equations.g # Value of all entries in diagonal of A31 (diagonal matrix).
    A13 = zeros(FT, polydeg, polydeg, nz, NumG)
    A23 = zeros(FT, polydeg, polydeg, nz, NumG)
    A32 = zeros(FT, polydeg, polydeg, nz, NumG)
    B1m_row1_23 = zeros(FT, 2, nz, NumG) # Col entries 2 and 3 of row 1 of first (upper) B block connecting each element with the one on the left
    B2m_row1_23 = zeros(FT, 2, nz, NumG) # Col entries 2 and 3 of row 1 of second (middle) B block connecting each element with the one on the left
    B3m_row1_23 = zeros(FT, 2, nz, NumG) # Col entries 2 and 3 of row 1 of last (lower) B block connecting each element with the one on the left
    B1_col2 = zeros(FT, polydeg, nz, NumG) # First (upper) B block connecting each elem with itself: polydeg (column 2)
    B2_col2 = zeros(FT, polydeg, nz, NumG) # Second (middle) B block connecting each elem with itself: polydeg (column 2)
    B3_col3 = zeros(FT, polydeg, nz, NumG) # Third (lower) B block connecting each elem with itself: polydeg (column 3)
    C2_row2 = zeros(FT, polydeg, nz, NumG) 
    C3_rows13 = zeros(FT, 2, polydeg, nz, NumG) 
    C2p_col1 = zeros(FT, 3, nz, NumG) 
    C3p_col1 = zeros(FT, 3, nz, NumG) 
    SA = zeros(FT, polydeg, polydeg, nz, NumG)
    SchurD = zeros(FT,3,3,nz,NumG)
    SchurL = zeros(FT,3,2,nz-1,NumG)
    SchurU =  zeros(FT,3,2,nz-1,NumG)
    Schurm =  zeros(FT,nz,NumG)
    Schurd1inv = zeros(FT,nz,NumG)
    Schurd1vw = zeros(FT,2,nz,NumG)
    Schuru1 = zeros(FT,2,nz-1,NumG)
    Schurl1 = zeros(FT,2,nz-1,NumG)
    SchurLhat = zeros(FT,2,2,nz-1,NumG)
    SchurDinv = zeros(FT,2,2,nz,NumG)
    SchurUhat = zeros(FT,2,2,nz-1,NumG)

    rs = zeros(FT, 3 * nz, NumG)

    polydeg = M - 1
    zCol = compute_vertical_step_size(NumG, polydeg, semi, nx, nz,
                                      typeof(semi.mesh))
    wPos, ThPos = compute_variable_position(equations)
    
    jacobian_aux = [JacobianAuxAndres(M, FT) for _ in 1:Threads.nthreads()]
    rotation_matrix = pre_compute_rotation_matrix(semi, equations)
    ID_map, iz_map = build_vertical_mapping(semi, nx, nz, equations)
    u_vert = zeros(FT, M, nz, NumG, nvariables(semi))
    return JacobianCacheAndres{nz, M, FT,
                         typeof(A12_11),
                         typeof(B1m_row1_23),
                         typeof(A13), typeof(zCol), typeof(jacobian_aux),
                         typeof(rotation_matrix), typeof(ID_map), typeof(iz_map),
                         typeof(u_vert)}(M, nz,
                                         nx,
                                         fac,
                                         A12_11,
                                         A22_diag,
                                         A33_diag,
                                         FacGrav,
                                         A13,
                                         A23,
                                         A32,
                                         B1m_row1_23,
                                         B2m_row1_23,
                                         B3m_row1_23,
                                         B1_col2,
                                         B2_col2,
                                         B3_col3,
                                         C2_row2,
                                         C3_rows13,
                                         C2p_col1,
                                         C3p_col1,
                                         SA,
                                        SchurD,
                                        SchurL,
                                        SchurU,
                                        Schurm,
                                        Schurd1inv,
                                        Schurd1vw,
                                        Schuru1,
                                        Schurl1,  
                                        SchurLhat,
                                        SchurDinv,
                                        SchurUhat,
                                         rs,
                                         NumG,
                                         zCol,
                                         cS,
                                         jacobian_aux,
                                         wPos,
                                         ThPos,
                                         rotation_matrix, ID_map, iz_map, u_vert)
end

@inline nvelem(::JacobianCacheAndres{NZ, M}) where {NZ, M} = NZ
@inline nvnodes(::JacobianCacheAndres{NZ, M}) where {NZ, M} = M

include("block_banded_matrix_v3.jl")

struct JacobianAuxAndres{R2Type, ThType, SALType}
    # SchurBoundary
    r2::R2Type
    r3::R2Type
    Th::ThType
    dpdRhoTh::ThType
    SAL::SALType
end

function JacobianAuxAndres(M, FT)
    r2 = zeros(FT, M - 1, 2)  # Auxiliary array to store polydeg*2 values (usually for rho_theta rhs)
    r3 = zeros(FT, M - 1, 2)  # Auxiliary array to store polydeg*2 values (usually for rho_w rhs)
    Th = zeros(FT, M)         # Auxiliary array to store polydeg+1 values (usually for rho_theta DOFs)
    dpdRhoTh = zeros(FT, M)   # Auxiliary array to store polydeg+1 values (usually for dp/d(rho_theta) DOFs)
    SAL = zeros(FT, M-1, M-1) # Auxiliary matrix of size polydeg x polydeg (for local Schur computations)
    return JacobianAuxAndres{typeof(r2), typeof(Th),
                             typeof(SAL)}(r2, r3, Th, dpdRhoTh, SAL)
end

@muladd begin
    # This function fills cache_jacobian::JacobianCacheAndres and computes the small Schur complements (SA) of each diagonal block of the A matrix, which is the upper-left block of the Jacobian matrix that corresponds to each column of LGL nodes.
    # Attention:
    # * This assumes an affine mapping in the vertical direction 
    # * Here, we store D in Schur
    function FillJacDGVertKernel!(cache_jacobian::JacobianCacheAndres, u, semi, equations)
        (; A13, A23, A32, A12_11, A22_diag, A33_diag, B1m_row1_23, B2m_row1_23, B3m_row1_23,  B1_col2, B2_col2, B3_col3, C2_row2, C3_rows13, C2p_col1, C3p_col1, SA, SchurD, SchurL, SchurU, M, jacobian_aux, fac, FacGrav, wPos, ThPos) = cache_jacobian

        (; nodes, derivative_hat, weights) = semi.solver.basis
        (; NumG, dz, cS) = cache_jacobian
        kappa = equations.R / equations.c_p
        invcS = 1 / cS
        kexp = kappa / (eltype(u)(1) - kappa)
        kfac = eltype(u)(1) / (eltype(u)(1) - kappa) * equations.R
        pre_fac = kfac * (equations.R / equations.p_0)^kexp
        DWS = derivative_hat
        DW = derivative_hat
        wB = weights[1]
        invwB = inv(wB)
        RhoPos = 1
        invfac = inv(fac)
        nz = nvelem(cache_jacobian)
        @threaded for ID in 1:NumG
            aux = jacobian_aux[Threads.threadid()]
            Th = aux.Th
            dpdRhoTh = aux.dpdRhoTh
            SAL = aux.SAL
            @inbounds for iz in 1:nz
                inv2dz = eltype(u)(2) / dz[iz, ID] # Affine mapping scaling
                invdz = eltype(u)(1) / dz[iz, ID]
                @inbounds for i in 1:M
                    Th[i] = u[i, iz, ID, ThPos] / u[i, iz, ID, RhoPos]
                    #dpdRhoTh[i] = pre_fac * exp(kexp * log(u[i, iz, ID, ThPos]))
                    dpdRhoTh[i] = pre_fac * u[i, iz, ID, ThPos]^kexp
                end
                
                # A block
                #########

                # Fill block A13 with internal contribution to d(rho_t)/d(rhow) 
                # (just derivative matrix with scaling)
                @inbounds for j in 1:(M - 1)
                    @inbounds for i in 1:(M - 1)
                        A13[i, j, iz, ID] = inv2dz * DWS[i, j]
                    end
                end
                # Zero first entry for elements iz=2, ..., nz
                # (Due to surface flux.. Modify for other surface fluxes?)
                if iz ≠ 1
                    A13[1, 1, iz, ID] = eltype(u)(0)
                end

                # Now fill d(rhotheta_t)/d(rhow) and d(rhow_t)/d(rhow)
                @inbounds for j in 1:(M - 1)
                    @inbounds for i in 1:(M - 1)
                        A23[i, j, iz, ID] = inv2dz * DWS[i, j] * Th[j]
                        A32[i, j, iz, ID] = inv2dz * DWS[i, j] * dpdRhoTh[j]
                    end
                end

                # Zero first entries for elements iz=2, ..., nz
                # (Due to surface flux.. Modify for other surface fluxes?)
                if iz ≠ 1
                    A23[1, 1, iz, ID] = eltype(u)(0)
                    A32[1, 1, iz, ID] = eltype(u)(0)
                else
                    A32[1, 1, iz, ID] *= -one(eltype(u))
                end

                # Now fill the additional entries that we don't have in Oswald's A
                # A12_11 = d(rho_t)/d(rhotheta) for node on "left" end of element
                # A22_11 = d(rhotheta_t)/d(rhotheta) for node on "left" end of element
                # A33_11 = d(rhow_t)/d(rhow) for node on "left" end of element (could be done in preprocessing)
                if iz == 1
                    # A12_11[iz, ID] = eltype(u)(0) # Already done in preprocessing
                    # A22_diag[1, iz, ID] = fac # Already done in preprocessing
                    A33_diag[1, iz, ID] = fac + inv2dz * cS * invwB
                else
                    A12_11[iz, ID] = dpdRhoTh[1] * invcS * invdz * invwB
                    A22_diag[1, iz, ID] = fac + Th[1] * dpdRhoTh[1] * invcS * invdz * invwB # TODO: We could also save 1/diag
                    A33_diag[1, iz, ID] = fac + cS * invdz * invwB
                end

                # Compute Schur complement for blocks of A
                @inbounds for i in 1:(M - 1)
                    @inbounds for j in 1:(M - 1)
                        SAL[i, j] = eltype(u)(0)

                        @inbounds for k in 1:(M - 1)
                            SAL[i, j] += A32[i, k, iz, ID] * A23[k, j, iz, ID] / A22_diag[k, iz, ID]
                        end
                        if i == j
                            SA[i, j, iz, ID] = A33_diag[i, iz, ID] -
                                               FacGrav * invfac * A13[i, j, iz, ID] -
                                               SAL[i, j]
                        else
                            SA[i, j, iz, ID] = -FacGrav * invfac * A13[i, j, iz, ID] -
                                               SAL[i, j]
                        end
                    end
                end
                # Add off-diagonal entry's correction
                @inbounds for j in 1:(M - 1)
                    SA[1, j, iz, ID] += (FacGrav * invfac * A23[1, j, iz, ID] * 
                                         A12_11[iz, ID] / A22_diag[1, iz, ID])
                end
                LUFull!(iz, ID, SA, Val(nvnodes(cache_jacobian) - 1))

                @inbounds for i in 1:(M-1)
                    B1_col2[i, iz, ID] = eltype(u)(2) * DW[i, M] / dz[iz, ID]
                    B2_col2[i, iz, ID] = inv2dz * DWS[i, M] * Th[M]
                    B3_col3[i, iz, ID] = inv2dz * DWS[i, M] * dpdRhoTh[M]
                end
                
                # B block
                #########

                # Surface flux entries (LMARS) for B matrix (connecting each element with the one on the "left")
                if iz > 1
                    # Density
                    B1m_row1_23[1, iz, ID] = -eltype(u)(1) / (wB * dz[iz, ID])
                    #dpdRhoThM = pre_fac * exp(kexp * log(u[M, iz - 1, ID, ThPos]))
                    dpdRhoThM = pre_fac * u[M, iz - 1, ID, ThPos]^kexp
                    B1m_row1_23[2, iz, ID] = -dpdRhoThM * invcS * invwB / dz[iz, ID]
                    # Potential temperature
                    ThM = u[M, iz - 1, ID, ThPos] / u[M, iz - 1, ID, RhoPos]
                    B2m_row1_23[1, iz, ID] = -ThM * invdz * invwB
                    B2m_row1_23[2, iz, ID] = -Th[1] * dpdRhoThM * invcS * invdz * invwB
                    # Momentum
                    B3m_row1_23[1, iz, ID] = -cS * invdz * invwB
                    B3m_row1_23[2, iz, ID] = -dpdRhoThM * invdz * invwB
                end

                # C block
                #########
                
                @inbounds for i in 1:(M-1)
                    # d(rhow_t)/d(rhotheta) for right DOF of each element wrt nodes 1:polydeg
                    C2_row2[i, iz, ID] = inv2dz * DWS[M, i] * dpdRhoTh[i] 

                    # d(rho_t)/d(rhow) for right DOF of each element wrt nodes 1:polydeg
                    C3_rows13[1, i, iz, ID] = inv2dz * DWS[M, i]
                    # d(rhotheta_t)/d(rhow) for right DOF of each element wrt nodes 1:polydeg
                    C3_rows13[2, i, iz, ID] = inv2dz * DWS[M, i] * Th[i]
                end

                if iz > 1
                    dpdRhoThP = pre_fac * u[1, iz, ID, ThPos]^kexp
                    ThM = u[M, iz - 1, ID, ThPos] / u[M, iz - 1, ID, RhoPos]
                
                    C2p_col1[1, iz - 1, ID] = -dpdRhoThP * invcS * invwB / dz[iz - 1, ID]
                    C2p_col1[2, iz - 1, ID] = dpdRhoTh[1] * invwB / dz[iz - 1, ID]
                    C2p_col1[3, iz - 1, ID] = -ThM * dpdRhoTh[1] * invcS * invwB / dz[iz - 1, ID]

                    C3p_col1[1, iz - 1, ID] = invwB / dz[iz - 1, ID]
                    C3p_col1[2, iz - 1, ID] = -cS * invwB / dz[iz - 1, ID]
                    C3p_col1[3, iz - 1, ID] = Th[1] * invwB / dz[iz - 1, ID]
                end

                # D block
                #########

                D11 = fac
                D21 = FacGrav
                D31 = eltype(u)(0)

                if iz < nz
                    ThM = u[M, iz, ID, ThPos] / u[M, iz, ID, RhoPos]
                    dpdRhoThM = pre_fac * u[M, iz, ID, ThPos]^kexp
                    D12 = eltype(u)(0)
                    D22 = fac + cS / dz[iz, ID] * invwB
                    D32 = eltype(u)(0)
                    D13 = dpdRhoTh[M] * invcS * invdz * invwB
                    D23 = eltype(u)(0)
                    D33 = fac + ThM * dpdRhoThM * invcS * invwB / dz[iz, ID]
                else
                    D12 = eltype(u)(2) * DW[M, M] / dz[iz, ID]
                    D22 = fac + inv2dz * cS * invwB
                    D32 = inv2dz * DWS[M, M] * Th[M]
                    D13 = eltype(u)(0)
                    D23 = -inv2dz * DWS[M, M] * dpdRhoTh[M]
                    D33 = fac
                end

                SchurD[1, 1, iz, ID] = D11
                SchurD[2, 1, iz, ID] = D21
                SchurD[3, 1, iz, ID] = D31
                SchurD[1, 2, iz, ID] = D12
                SchurD[2, 2, iz, ID] = D22
                SchurD[3, 2, iz, ID] = D32
                SchurD[1, 3, iz, ID] = D13
                SchurD[2, 3, iz, ID] = D23
                SchurD[3, 3, iz, ID] = D33

                # Reset lower and upper blocks
                if iz < nz
                    f0 = eltype(u)(0)
                    SchurL[1, 1, iz, ID] = f0
                    SchurL[2, 1, iz, ID] = f0
                    SchurL[3, 1, iz, ID] = f0
                    SchurL[1, 2, iz, ID] = f0
                    SchurL[2, 2, iz, ID] = f0
                    SchurL[3, 2, iz, ID] = f0
                    SchurU[1, 1, iz, ID] = f0
                    SchurU[2, 1, iz, ID] = f0
                    SchurU[3, 1, iz, ID] = f0
                    SchurU[1, 2, iz, ID] = f0
                    SchurU[2, 2, iz, ID] = f0
                    SchurU[3, 2, iz, ID] = f0
                end
            end # end ID
        end # end iz
    end
end

"""
    SchurBoundaryKernel!(cache::JacobianCacheAndres, equations)

Compute  Schur ← D − C·A⁻¹·B in-place for each vertical column.
`Schur` must be pre-loaded with D values before this call.

### SA Schur complement  (eliminates ρw = variable 3)
  For each diagonal block of A:
  SA = A33 − A31·A11⁻¹·A13 − A32·A22⁻¹·A23  (pre-factored, Doolittle LU)
  SA rhs for B column [b1; b2; b3]:
    rhs[i] = b3·δᵢ₁ − FacGrav·invfac·b1·δᵢ₁ − Σⱼ A32[i,j]·b2[j]/A22_diag[j]
            + FacGrav·invfac·A12_11·b2[1]/A22_diag[1]·δᵢ₁      ← A12 cross-term

### Off-diagonal blocks  (iz > 1 only)
  Below-diagonal B_m (B1m/B2m/B3m_row1_23, only interior row 1):
    ρ/ρθ/ρw interior node-1 of iz driven by edge cols 2,3 of iz-1
    → sub-diagonal  S[iz, iz-1] = Schur.L[iz-1]   (two full SA solves)
  Above-diagonal C_m (C2p_col1, C3p_col1, only interior col 1):
    C2p_col1[j]: C[edge eq j of iz, ρθ-interior node-1 of iz+1]
    C3p_col1[j]: C[edge eq j of iz, ρw-interior node-1 of iz+1]
    → super-diagonal  S[iz-1, iz] = Schur.U[iz-1]  (O(1) with saved row-1 values)

### Key optimisations
1. **Batched 2-column SA solve** — cols 2 & 3 (diagonal) share one forward/
   backward pass; km=1 & km=2 (sub-diagonal) share another.  SA is read
   **twice** per element instead of four times.
2. **Column-major A23 / A32 loops** — outer j, inner i; `A23[:,j,iz,ID]`
   and `A32[:,j,iz,ID]` are contiguous columns → stride-1 `@simd` loads.
3. **All O(polydeg²) ops process both columns at once** — A23, A32,
   C3_rows13, C2_row2, A13 each traversed once per element.
4. **Column views** — `@view r3[:, 1]` / `r3[:, 2]` at the top of each iz
   iteration give the vectoriser unambiguous non-aliasing 1-D information.

### Auxiliary struct
`jacobian_aux[tid]` needs only **two** 2-D scratch arrays:
    r2 = zeros(FT, M - 1, 2)   # ρθ interior; columns 1/2 = col-3/col-2 (or km1/km2)
    r3 = zeros(FT, M - 1, 2)   # ρw interior; same convention

Column-major layout makes `r3[:, 1]` and `r3[:, 2]` co-located in memory
(polydeg * sizeof(FT) apart), which is marginally better than two separate
1-D arrays for cache warm-up on small polydeg.
"""
function SchurBoundaryKernel!(cache::JacobianCacheAndres, equations)
    (; B1m_row1_23, B2m_row1_23, B3m_row1_23,
       B1_col2, B2_col2, B3_col3,
       C2_row2, C3_rows13, C2p_col1, C3p_col1,
       A13, A23, A32, A12_11, A22_diag,
       SA, SchurD, SchurL, SchurU, M, NumG, jacobian_aux, fac, FacGrav) = cache
 
    FT      = eltype(SA)
    invfac  = inv(fac)
    nz      = cache.nz
    polydeg = M - 1
 
    @threaded for ID in 1:NumG
        aux = jacobian_aux[Threads.threadid()]
        r2  = aux.r2    # (polydeg, 2)
        r3  = aux.r3    # (polydeg, 2)
 
        @inbounds for iz in 1:nz
 
            r3a = @view r3[:, 1];   r3b = @view r3[:, 2]
            r2a = @view r2[:, 1];   r2b = @view r2[:, 2]
 
            a12 = A12_11[iz, ID]
            fgi = FacGrav * invfac
 
            # ── Step 1: A22⁻¹·B2_col2 → r2a ────────────────────────────
            @simd for i in 1:polydeg
                r2a[i] = B2_col2[i, iz, ID] / A22_diag[i, iz, ID]
            end
 
            # ── Step 2: build both SA rhs ────────────────────────────────
            # r3a ← B3_col3  (col-3 / ρθ-edge batch)
            # r3b ← −(FacGrav·invfac·B1 + A32·r2a)  (col-2 / ρw-edge batch)
            #   outer j: NO @simd — r3b accumulates across j
            #   inner i: @simd — A32[:,j] contiguous, independent across i
            @simd for i in 1:polydeg
                r3a[i] = B3_col3[i, iz, ID]
                r3b[i] = -fgi * B1_col2[i, iz, ID]
            end
            for j in 1:polydeg
                r2aj = r2a[j]
                @simd for i in 1:polydeg
                    r3b[i] -= A32[i, j, iz, ID] * r2aj
                end
            end
            r3b[1] += fgi * a12 * r2a[1]
 
            # ── Step 3: batched 2-column SA solve [r3a | r3b] ───────────
            #   outer k: NO @simd;  inner i: @simd
            for k in 1:polydeg-1
                r3ak = r3a[k];   r3bk = r3b[k]
                @simd for i in k+1:polydeg
                    sa_ik  = SA[i, k, iz, ID]
                    r3a[i] -= sa_ik * r3ak
                    r3b[i] -= sa_ik * r3bk
                end
            end
            for k in polydeg:-1:1
                sa_inv = inv(SA[k, k, iz, ID])
                r3a[k] *= sa_inv;   r3b[k] *= sa_inv
                r3ak = r3a[k];   r3bk = r3b[k]
                @simd for i in 1:k-1
                    sa_ik  = SA[i, k, iz, ID]
                    r3a[i] -= sa_ik * r3ak
                    r3b[i] -= sa_ik * r3bk
                end
            end
            r3A1 = r3a[1];   r3B1 = r3b[1]   # saved for super-diagonal
 
            # ── Step 4: batched A23 → [r2a, r2b] ────────────────────────
            #   init: @simd;  outer j: NO @simd;  inner i: @simd;  divide: @simd
            @simd for i in 1:polydeg
                r2a[i] = zero(FT)
                r2b[i] = B2_col2[i, iz, ID]
            end
            for j in 1:polydeg
                r3aj = r3a[j];   r3bj = r3b[j]
                @simd for i in 1:polydeg
                    aij    = A23[i, j, iz, ID]
                    r2a[i] -= aij * r3aj
                    r2b[i] -= aij * r3bj
                end
            end
            @simd for i in 1:polydeg
                ia22   = inv(A22_diag[i, iz, ID])
                r2a[i] *= ia22;   r2b[i] *= ia22
            end
            r2A1 = r2a[1];   r2B1 = r2b[1]   # saved for super-diagonal
 
            # ── Step 5: batched ρ node-1 recovery (scalar reductions) ───
            r1A = -a12 * r2a[1]
            r1B =  B1_col2[1, iz, ID] - a12 * r2b[1]
            @simd for j in 1:polydeg
                a13j = A13[1, j, iz, ID]
                r1A -= a13j * r3a[j]
                r1B -= a13j * r3b[j]
            end
            r1A *= invfac;   r1B *= invfac
 
            # ── Step 6: batched C_self correction → [s3, s2] ────────────
            # s accumulators use −=, so s = −(C·A⁻¹·B).
            # Write as D + s (not D − s) to get D − C·A⁻¹·B.
            s3_1 = zero(FT);  s3_2 = zero(FT);  s3_3 = zero(FT)
            s2_1 = zero(FT);  s2_2 = zero(FT);  s2_3 = zero(FT)
            @simd for j in 1:polydeg
                c1j = C3_rows13[1, j, iz, ID]
                c3j = C3_rows13[2, j, iz, ID]
                c2j = C2_row2[j,    iz, ID]
                s3_1 -= c1j * r3a[j];   s3_3 -= c3j * r3a[j];   s3_2 -= c2j * r2a[j]
                s2_1 -= c1j * r3b[j];   s2_3 -= c3j * r3b[j];   s2_2 -= c2j * r2b[j]
            end
 
            # ── Write D[iz] ─────────────────────────────────────────
            # D + s  (s already carries the minus sign from accumulation).
            SchurD[1, 2, iz, ID] += s2_1
            SchurD[2, 2, iz, ID] += s2_2
            SchurD[3, 2, iz, ID] += s2_3
            SchurD[1, 3, iz, ID] += s3_1
            SchurD[2, 3, iz, ID] += s3_2
            SchurD[3, 3, iz, ID] += s3_3
 
            if iz > 1
 
                # ── Super-diagonal U[iz-1]: O(1) ───────────────────
                # Direct subtraction of C_p · (node-1 diagonal solutions).
                SchurU[1, 1, iz - 1, ID] -= C2p_col1[1,iz-1,ID]*r2B1 + C3p_col1[1,iz-1,ID]*r3B1
                SchurU[2, 1, iz - 1, ID] -= C2p_col1[2,iz-1,ID]*r2B1 + C3p_col1[2,iz-1,ID]*r3B1
                SchurU[3, 1, iz - 1, ID] -= C2p_col1[3,iz-1,ID]*r2B1 + C3p_col1[3,iz-1,ID]*r3B1
                SchurU[1, 2, iz - 1, ID] -= C2p_col1[1,iz-1,ID]*r2A1 + C3p_col1[1,iz-1,ID]*r3A1
                SchurU[2, 2, iz - 1, ID] -= C2p_col1[2,iz-1,ID]*r2A1 + C3p_col1[2,iz-1,ID]*r3A1
                SchurU[3, 2, iz - 1, ID] -= C2p_col1[3,iz-1,ID]*r2A1 + C3p_col1[3,iz-1,ID]*r3A1
 
                # ── Sub-diagonal B_m scalars ──────────────────────────────
                b1m1 = B1m_row1_23[1, iz, ID];   b1m2 = B1m_row1_23[2, iz, ID]
                b2m1 = B2m_row1_23[1, iz, ID];   b2m2 = B2m_row1_23[2, iz, ID]
                b3m1 = B3m_row1_23[1, iz, ID];   b3m2 = B3m_row1_23[2, iz, ID]
                ia22_1  = inv(A22_diag[1, iz, ID])
                b2m1_sc = b2m1 * ia22_1;   b2m2_sc = b2m2 * ia22_1
 
                # ── Sub-diagonal SA rhs for km=1,2 ───────────────────────
                @simd for i in 1:polydeg
                    a32i1  = A32[i, 1, iz, ID]
                    r3a[i] = -b2m1_sc * a32i1
                    r3b[i] = -b2m2_sc * a32i1
                end
                r3a[1] += b3m1 - fgi * b1m1 + fgi * a12 * b2m1_sc
                r3b[1] += b3m2 - fgi * b1m2 + fgi * a12 * b2m2_sc
 
                # ── Batched SA solve for km=1,2 ───────────────────────────
                for k in 1:polydeg-1
                    r3ak = r3a[k];   r3bk = r3b[k]
                    @simd for i in k+1:polydeg
                        sa_ik  = SA[i, k, iz, ID]
                        r3a[i] -= sa_ik * r3ak
                        r3b[i] -= sa_ik * r3bk
                    end
                end
                for k in polydeg:-1:1
                    sa_inv = inv(SA[k, k, iz, ID])
                    r3a[k] *= sa_inv;   r3b[k] *= sa_inv
                    r3ak = r3a[k];   r3bk = r3b[k]
                    @simd for i in 1:k-1
                        sa_ik  = SA[i, k, iz, ID]
                        r3a[i] -= sa_ik * r3ak
                        r3b[i] -= sa_ik * r3bk
                    end
                end
 
                # ── A23 recovery for km=1,2 ───────────────────────────────
                r2a[1] = b2m1;   r2b[1] = b2m2
                @simd for i in 2:polydeg
                    r2a[i] = zero(FT);   r2b[i] = zero(FT)
                end
                for j in 1:polydeg
                    r3aj = r3a[j];   r3bj = r3b[j]
                    @simd for i in 1:polydeg
                        aij    = A23[i, j, iz, ID]
                        r2a[i] -= aij * r3aj
                        r2b[i] -= aij * r3bj
                    end
                end
                @simd for i in 1:polydeg
                    ia22   = inv(A22_diag[i, iz, ID])
                    r2a[i] *= ia22;   r2b[i] *= ia22
                end
 
                # ── C_p correction to D[iz-1] ─────────────────────────────
                # S[iz-1,iz-1] also receives −C_p[iz-1]·A⁻¹[iz]·B_m[iz,iz-1].
                # A⁻¹[iz]·B_m is already solved above; only node 1 enters
                # (C_p is sparse at col 1).  Direct subtraction — no sign issue.  ← FIX 3
                x2a1 = r2a[1];   x3a1 = r3a[1]   # km=1 (col 2 / ρw-edge)
                x2b1 = r2b[1];   x3b1 = r3b[1]   # km=2 (col 3 / ρθ-edge)
 
                # ── ρ node-1 recovery for km=1,2 (scalar reductions) ─────
                r1m1 = b1m1 - a12 * r2a[1]
                r1m2 = b1m2 - a12 * r2b[1]
                @simd for j in 1:polydeg
                    a13j = A13[1, j, iz, ID]
                    r1m1 -= a13j * r3a[j]
                    r1m2 -= a13j * r3b[j]
                end
                r1m1 *= invfac;   r1m2 *= invfac

                SchurD[1, 2, iz - 1, ID] -= C2p_col1[1,iz-1,ID]*x2a1 + C3p_col1[1,iz-1,ID]*x3a1
                SchurD[2, 2, iz - 1, ID] -= C2p_col1[2,iz-1,ID]*x2a1 + C3p_col1[2,iz-1,ID]*x3a1
                SchurD[3, 2, iz - 1, ID] -= C2p_col1[3,iz-1,ID]*x2a1 + C3p_col1[3,iz-1,ID]*x3a1
                SchurD[1, 3, iz - 1, ID] -= C2p_col1[1,iz-1,ID]*x2b1 + C3p_col1[1,iz-1,ID]*x3b1
                SchurD[2, 3, iz - 1, ID] -= C2p_col1[2,iz-1,ID]*x2b1 + C3p_col1[2,iz-1,ID]*x3b1
                SchurD[3, 3, iz - 1, ID] -= C2p_col1[3,iz-1,ID]*x2b1 + C3p_col1[3,iz-1,ID]*x3b1
 
                # ── C_self correction for km=1,2 (scalar reductions) ─────
                # l accumulators use −=, so l = −(C·A⁻¹·B_m).
                # Write as L + l (not L − l) to get 0 − C·A⁻¹·B_m.  ← FIX 2
                l11 = zero(FT);  l21 = zero(FT);  l31 = zero(FT)
                l12 = zero(FT);  l22 = zero(FT);  l32 = zero(FT)
                @simd for j in 1:polydeg
                    c1j = C3_rows13[1, j, iz, ID]
                    c3j = C3_rows13[2, j, iz, ID]
                    c2j = C2_row2[j,    iz, ID]
                    l11 -= c1j * r3a[j];   l12 -= c1j * r3b[j]
                    l31 -= c3j * r3a[j];   l32 -= c3j * r3b[j]
                    l21 -= c2j * r2a[j];   l22 -= c2j * r2b[j]
                end
 
                # ── Write L[iz-1] ───────────────────────────────────
                SchurL[1, 1, iz - 1,  ID] += l11
                SchurL[2, 1, iz - 1,  ID] += l21
                SchurL[3, 1, iz - 1,  ID] += l31
                SchurL[1, 2, iz - 1,  ID] += l12
                SchurL[2, 2, iz - 1,  ID] += l22
                SchurL[3, 2, iz - 1,  ID] += l32
 
            end  # iz > 1
 
        end  # iz
    end  # ID
end


# Perform LU factorization of the Schur complement. Since this matrix has
# 4 lower bands and 5 upper bands, a banded-matrix format is not ideal for this matrix with l=4, u=5. 
# We can treat it as a block-banded matrix with 3x3 blocks.
function luBandkernel!(cache_jacobian::JacobianCacheAndres)
    block_lu!(cache_jacobian) # @threaded inside   
end

"""
    ldivVerticalFKernel!(cache::JacobianCacheAndres, b, equations)

Compute  rs = b_D − C·A⁻¹·b_A  for every column, storing the result in
`cache.rs` (size 3·nz × NumG, edge ordering ρ_R / ρw_R / ρθ_R per element).

### Derivation of the SA rhs
The interior block system is:

    [fac·I   A12   A13] [x1]   [b1]          (ρ,   nodes 1..M-1)
    [  0     A22   A23] [x2] = [b2]          (ρθ,  nodes 1..M-1)
    [FacGrav·I A32 A33] [x3]   [b3]          (ρw,  nodes 1..M-1)

Eliminating x2 from row 3 via row 2, then x1 via row 1:

    S32 = (A32 − FacGrav·invfac·A12) · A22⁻¹

    SA·x3 = b3 − FacGrav·invfac·b1 − S32·b2

Since A12 only has A12_11 ≠ 0 at [1,1], S32·b2 differs from A32·(b2/A22_diag)
only at row 1:

    rhs3[i] = b3[i] − FacGrav·invfac·b1[i] − Σⱼ A32[i,j]·b2[j]/A22_diag[j]
    rhs3[1] += FacGrav·invfac·A12_11·b2[1]/A22_diag[1]            (A12 correction)

After the SA solve (x3 ← SA⁻¹·rhs3):

    x2[i] = (b2[i] − Σⱼ A23[i,j]·x3[j]) / A22_diag[i]
    x1[1] = (b1[1] − A12_11·x2[1] − Σⱼ A13[1,j]·x3[j]) / fac    (scalar only)

### C correction
Self (element iz → right edge of iz):
    rs[sh+1] -= Σⱼ C3_rows13[1,j,iz]·x3[j]    (ρ_R   ← ρw interior)
    rs[sh+2] -= Σⱼ C2_row2[j,iz]·x2[j]         (ρw_R  ← ρθ interior)
              + FacGrav·x1[1]                    (ρw_R  ← gravity·ρ_interior[1])
    rs[sh+3] -= Σⱼ C3_rows13[2,j,iz]·x3[j]    (ρθ_R  ← ρw interior)

Above-diagonal (interior of iz → right edge of iz-1, node 1 only):
    rs[sh_prev+j] -= C3p_col1[j,iz-1]·x3[1] + C2p_col1[j,iz-1]·x2[1],  j=1,2,3

`C2p_col1[j,iz,ID]` = d(edge eq j of element iz)/d(ρθ at col 1 of element iz+1).
`C3p_col1[j,iz,ID]` = d(edge eq j of element iz)/d(ρw at col 1 of element iz+1).
Index is iz-1 (not iz) when correcting rs[iz-1] while the loop is at iz.

### Annotation conventions
@inbounds: one annotation on the `for iz` loop; propagates to all nested accesses.
@simd:     inner i-loops (elementwise/BLAS inner) and scalar reductions over j.
           NOT on outer j-loops (r2/r3 accumulate across j).
"""
function ldivVerticalFKernel!(cache_jacobian::JacobianCacheAndres, b, equations)
    (; M, A32, A13, A23, C3_rows13, C2_row2, C3p_col1, C2p_col1,
       A12_11, A22_diag, rs, SA, jacobian_aux, NumG, ThPos, wPos) = cache_jacobian

    FT      = eltype(SA)
    invfac  = inv(cache_jacobian.fac)
    FacGrav = equations.g
    nz      = nvelem(cache_jacobian)
    polydeg = M - 1
    RhoPos  = 1

    @threaded for ID in 1:NumG
        aux = jacobian_aux[Threads.threadid()]
        # Views taken once per thread, outside the iz hot loop.
        # @view on a column of a Matrix produces a contiguous SubArray (no heap alloc).
        r2 = @view aux.r2[:, 1]    # ρθ interior scratch
        r3 = @view aux.r3[:, 1]    # ρw interior scratch

        @inbounds for iz in 1:nz
            sh  = (iz - 1) * 3
            fgi = FacGrav * invfac

            # ── Init rs from right-edge b values (D-block seed) ─────────
            rs[sh + 1, ID] = b[M, iz, ID, RhoPos]
            rs[sh + 2, ID] = b[M, iz, ID, wPos]
            rs[sh + 3, ID] = b[M, iz, ID, ThPos]

            # ── Load r2 ← b_ρθ; r3 ← b_ρw − FacGrav·invfac·b_ρ ────────
            # Folds the FacGrav subtraction into the initialisation pass.
            # b_ρ read directly; no r1 scratch needed.
            @simd for i in 1:polydeg
                r2[i] = b[i, iz, ID, ThPos]
                r3[i] = b[i, iz, ID, wPos] - fgi * b[i, iz, ID, RhoPos]
            end

            # ── SA rhs: r3 -= A32·(b_ρθ / A22_diag) ────────────────────
            # Column-major: j outer (NO @simd — r3 accumulates over j),
            #               @simd for i (A32[:,j] contiguous, independent)
            for j in 1:polydeg
                r2j_sc = r2[j] / A22_diag[j, iz, ID]
                @simd for i in 1:polydeg
                    r3[i] -= A32[i, j, iz, ID] * r2j_sc
                end
            end
            # A12 cross-term: corrects row 1 only
            r3[1] += fgi * A12_11[iz, ID] * r2[1] / A22_diag[1, iz, ID]

            # ── SA solve: r3 ← SA⁻¹·r3 = x3 ────────────────────────────
            ldivFull!(iz, ID, SA, r3,  Val(nvnodes(cache_jacobian) - 1))

            # ── Recover x2: r2 ← (b_ρθ − A23·x3) / A22_diag ───────────
            # j outer (NO @simd — r2 accumulates over j), @simd for i
            for j in 1:polydeg
                r3j = r3[j]
                @simd for i in 1:polydeg
                    r2[i] -= A23[i, j, iz, ID] * r3j
                end
            end
            @simd for i in 1:polydeg
                r2[i] /= A22_diag[i, iz, ID]
            end

            # ── Self C correction → rs[iz] ───────────────────────────────
            # Load into scalars first; write back once.
            rs1 = rs[sh + 1, ID]
            rs2 = rs[sh + 2, ID]
            rs3 = rs[sh + 3, ID]
            @simd for j in 1:polydeg              # scalar reductions → @simd ok
                c1j = C3_rows13[1, j, iz, ID]
                c3j = C3_rows13[2, j, iz, ID]
                c2j = C2_row2[j,    iz, ID]
                rs1 -= c1j * r3[j]
                rs3 -= c3j * r3[j]
                rs2 -= c2j * r2[j]
            end
            rs[sh + 1, ID] = rs1
            rs[sh + 2, ID] = rs2
            rs[sh + 3, ID] = rs3

            # ── C_p correction → rs[iz-1] ─────────────────────────────
            # C2p_col1[j, iz, ID] = d(edge eq j of element iz) / d(ρθ col-1 of iz+1)
            # C3p_col1[j, iz, ID] = d(edge eq j of element iz) / d(ρw col-1 of iz+1)
            # "Col 1" = left boundary of element iz+1 (LGL) = right boundary of iz.
            # When processing element iz, its col-1 interior (x2[1], x3[1]) feeds
            # into edge equations of iz-1 → cache index is iz-1, not iz.
            if iz > 1
                sh_prev = (iz - 2) * 3
                x3_1 = r3[1]
                x2_1 = r2[1]
                @simd for j in 1:3
                    rs[sh_prev + j, ID] -= C3p_col1[j, iz - 1, ID] * x3_1 +
                                           C2p_col1[j, iz - 1, ID] * x2_1
                end
            end

        end  # iz
    end  # ID
end

# Solve for LU factorized block-three-diagonal matrices
function ldivVerticalSKernel!(cache_jacobian::JacobianCacheAndres)
    (; rs) = cache_jacobian
    solve!(cache_jacobian, rs) # @threaded inside
end

"""
    ldivVerticalBKernel!(cache::JacobianCacheAndres, b, equations)

After the Schur block-tridiagonal solve has written the edge solutions into `rs`,
this function:
  1. Copies rs → b at right-boundary node M  (edge DOFs: ρ_R, ρw_R, ρθ_R)
  2. Subtracts B·y from the interior RHS stored in b
  3. Solves the interior block system for x1 (ρ), x2 (ρθ), x3 (ρw)
  4. Writes the interior solutions back into b

### B-matrix structure (new right-edge-only permutation)
Self (col 2 = ρw-edge, col 3 = ρθ-edge; col 1 = ρ-edge has no coupling):
  B1_col2[i]  : ρ interior ← ρw-edge of iz        (all i = 1..polydeg)
  B2_col2[i]  : ρθ interior ← ρw-edge of iz
  B3_col3[i]  : ρw interior ← ρθ-edge of iz

Below-diagonal (only interior row 1 nonzero; km=1 ≡ ρw-edge of iz-1, km=2 ≡ ρθ-edge):
  B1m_row1_23[km] : ρ  interior row 1 ← edge km of iz-1
  B2m_row1_23[km] : ρθ interior row 1 ← edge km of iz-1
  B3m_row1_23[km] : ρw interior row 1 ← edge km of iz-1

No above-diagonal B (left-edge DOFs do not exist in the new permutation).

### Interior solve (same SA structure as ldivVerticalFKernel!)
  SA·x3 = b_ρw_mod − FacGrav·invfac·b_ρ_mod − A32·A22⁻¹·b_ρθ_mod   (+ A12 correction)
  x2    = (b_ρθ_mod − A23·x3) / A22_diag
  x1    = (b_ρ_mod  − A12_11·x2[1]·e1 − A13·x3) / fac

b_ρ is modified in-place; no separate r1 scratch array is needed.

### Annotation conventions
  @inbounds : one annotation on `for iz`; propagates inward.
  @simd     : inner i-loops (independent iterations, stride-1 access).
              NOT on outer j-loops (b_ρ / r2 / r3 accumulate over j).
"""
function ldivVerticalBKernel!(cache_jacobian::JacobianCacheAndres, b, equations)
    # TODO: Check correctness
    (; A13, A23, A32, B1_col2, B2_col2, B3_col3,
       B1m_row1_23, B2m_row1_23, B3m_row1_23,
       A12_11, A22_diag, SA, rs, M, NumG, jacobian_aux, FacGrav, ThPos, wPos) = cache_jacobian

    FT      = eltype(SA)
    invfac  = inv(cache_jacobian.fac)
    nz      = nvelem(cache_jacobian)
    polydeg = M - 1
    RhoPos  = 1

    @threaded for ID in 1:NumG
        aux = jacobian_aux[Threads.threadid()]
        r2  = @view aux.r2[:, 1]    # ρθ interior scratch
        r3  = @view aux.r3[:, 1]    # ρw interior scratch

        @inbounds for iz in 1:nz
            sh  = (iz - 1) * 3
            fgi = FacGrav * invfac

            # ── Step 1: write edge solution rs → b at right boundary ─────
            b[M, iz, ID, RhoPos] = rs[sh + 1, ID]
            b[M, iz, ID, wPos]   = rs[sh + 2, ID]
            b[M, iz, ID, ThPos]  = rs[sh + 3, ID]

            # ── Steps 2-3: load interior b, subtract B_self·y_edge ───────
            # Load r2 (ρθ) and r3 (ρw) while subtracting self-edge B in one
            # pass.  b_ρ is modified in-place (no r1 scratch needed).
            # B_col1 = 0: ρ-edge DOF has no coupling to interior.
            rsh2 = rs[sh + 2, ID]   # ρw-edge of iz
            rsh3 = rs[sh + 3, ID]   # ρθ-edge of iz
            @simd for i in 1:polydeg
                r2[i]                  = b[i, iz, ID, ThPos] - B2_col2[i, iz, ID] * rsh2
                r3[i]                  = b[i, iz, ID, wPos]  - B3_col3[i, iz, ID] * rsh3
                b[i, iz, ID, RhoPos] -= B1_col2[i, iz, ID] * rsh2
            end

            # ── Subtract B_m·y_{iz-1}: below-diagonal, row 1 only ────────
            if iz > 1
                sh_prev = (iz - 2) * 3
                rsp2 = rs[sh_prev + 2, ID]   # ρw-edge of iz-1 (km=1)
                rsp3 = rs[sh_prev + 3, ID]   # ρθ-edge of iz-1 (km=2)
                r2[1] -= B2m_row1_23[1, iz, ID] * rsp2 + B2m_row1_23[2, iz, ID] * rsp3
                r3[1] -= B3m_row1_23[1, iz, ID] * rsp2 + B3m_row1_23[2, iz, ID] * rsp3
                b[1, iz, ID, RhoPos] -= B1m_row1_23[1, iz, ID] * rsp2 +
                                        B1m_row1_23[2, iz, ID] * rsp3
            end

            # ── Step 4: SA rhs → solve for ρw interior ───────────────────
            # r3 ← b_ρw_mod − FacGrav·invfac·b_ρ_mod − A32·(r2/A22_diag)
            # b[:,iz,RhoPos] holds the modified ρ interior RHS here.
            # @simd on gravity loop (elementwise); outer j NO @simd (r3 accumulates)
            @simd for i in 1:polydeg
                r3[i] -= fgi * b[i, iz, ID, RhoPos]
            end
            for j in 1:polydeg
                r2j_sc = r2[j] / A22_diag[j, iz, ID]
                @simd for i in 1:polydeg
                    r3[i] -= A32[i, j, iz, ID] * r2j_sc
                end
            end
            r3[1] += fgi * A12_11[iz, ID] * r2[1] / A22_diag[1, iz, ID]

            ldivFull!(iz, ID, SA, r3,  Val(nvnodes(cache_jacobian) - 1))   # r3 ← x3

            # ── Step 5: recover ρθ interior x2 → r2 ─────────────────────
            # x2 = (b_ρθ_mod − A23·x3) / A22_diag
            # outer j: NO @simd (r2 accumulates); inner @simd
            for j in 1:polydeg
                r3j = r3[j]
                @simd for i in 1:polydeg
                    r2[i] -= A23[i, j, iz, ID] * r3j
                end
            end
            @simd for i in 1:polydeg
                r2[i] /= A22_diag[i, iz, ID]
            end

            # ── Step 6: recover ρ interior x1, overwrite b_ρ in-place ───
            # x1 = (b_ρ_mod − A12_11·x2[1]·e1 − A13·x3) / fac
            # A12 correction at row 1 only; outer j: NO @simd (b_ρ accumulates)
            b[1, iz, ID, RhoPos] -= A12_11[iz, ID] * r2[1]
            for j in 1:polydeg
                r3j = r3[j]
                @simd for i in 1:polydeg
                    b[i, iz, ID, RhoPos] -= A13[i, j, iz, ID] * r3j
                end
            end
            @simd for i in 1:polydeg
                b[i, iz, ID, RhoPos] *= invfac
            end

            # ── Step 7: write ρθ and ρw solutions to b ───────────────────
            @simd for i in 1:polydeg
                b[i, iz, ID, ThPos] = r2[i]
                b[i, iz, ID, wPos]  = r3[i]
            end

        end  # iz
    end  # ID
end
