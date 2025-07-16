abstract type AbstractSolver end

mutable struct ElectronPhononSolver{MT, VT} <: AbstractSolver
    iter :: Int

    model :: MT
    η :: Float64
    occupation :: Float64

    qpts :: Kpoints{VT}

    # ωs and ks are a Vector so they can be nonuniform.
    # ωs_dense and ks_dense are a StepRangeLen so linear interpolation over them are efficient.

    # ωs :: StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    ωs :: Vector{Float64}
    ks :: Vector{Float64}
    ωs_dense :: StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    ks_dense :: StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    Σs :: Matrix{ComplexF64}
    As :: Matrix{Float64}
    Σs_dense :: Matrix{ComplexF64}
    As_dense :: Matrix{Float64}

    Σs_rest :: Vector{ComplexF64}

    spectral_sum :: Vector{Float64}
    spectral_occ :: Vector{Float64}

    # Precomputed integral weights
    # Kramers-Kronig transformation matrix for ωs (Re Σ = KK_map * Im Σ)
    KK_map :: Matrix{Float64}
    # Integration weight for f(k) = ∫d³k' g(k') / |k - k'|² becomes
    # f[ik] = ∑_ikq q_integral_matrix[ik, ikq] * g[ikq] where f[ik] and g[ik] are linear
    # spline coefficients
    q_integral_matrix :: Matrix{Float64}
end

function ElectronPhononSolver(model, ωs_, ks_, qpts; occupation, ωs_dense = nothing, ks_dense = nothing, η, Σs_rest = nothing)
    # Convert ranges to vectors
    ωs = Vector(ωs_)
    ks = Vector(ks_)

    # Default values for dense grids: same range, but 5 times denser
    if ωs_dense === nothing
        ωs_dense = range(extrema(ωs)..., step = minimum(diff(ωs)) / 5)
    end
    if ks_dense === nothing
        ks_dense = range(extrema(ks)..., step = minimum(diff(ks)) / 5)
    end

    Σs = zeros(ComplexF64, length(ωs), length(ks))
    As = zeros(length(ωs), length(ks))
    Σs_dense = zeros(ComplexF64, length(ωs_dense), length(ks_dense))
    As_dense = zeros(length(ωs_dense), length(ks_dense))
    spectral_sum = zeros(length(ks_dense))
    spectral_occ = zeros(length(ks_dense))

    if Σs_rest === nothing
        Σs_rest = zeros(ComplexF64, length(ks))
    end

    KK_map = get_Kramers_Kronig_map(ωs) ./ π

    # For efficient q integral, we use linear spline basis. (qpts is no longer used.)
    # This requires a uniform k mesh.
    # (For nonuniform k mesh, frohlich_compute_angle_factor must be updated.)
  #  if model isa FrohlichModel
  #      ks_tmp, _, q_integral_matrix = frohlich_compute_angle_factor(maximum(ks), length(ks), 1);
  #      q_integral_matrix .*= abs2(get_eph_g(1.0, model))
  #      @assert all(diff(ks) .≈ diff(ks)[1])
  #      @assert ks ≈ ks_tmp
  #  else
        q_integral_matrix = zeros(length(ks), length(ks))
  #  end

    ElectronPhononSolver(0, model, η, occupation, qpts, ωs, ks, ωs_dense, ks_dense,
        Σs, As, Σs_dense, As_dense, Σs_rest, spectral_sum, spectral_occ,
        KK_map, q_integral_matrix)
end

function Base.show(io :: IO, S :: ElectronPhononSolver)
    print(io, "ElectronPhononSolver(", S.model, ", η = ", S.η, ")\n")
    print(io, "qpts     : ", length(S.qpts), " points, ", "\n")
    print(io, "ωs       : ", length(S.ωs), " points, ", S.ωs, "\n")
    print(io, "ks       : ", length(S.ks), " points, ", S.ks, "\n")
    print(io, "ωs_dense : ", length(S.ωs_dense), " points, ", S.ωs_dense, "\n")
    print(io, "ks_dense : ", length(S.ks_dense), " points, ", S.ks_dense)
end