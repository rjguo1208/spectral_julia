function get_Σ_itp(S :: PhononDefectSolver, η = S.η)
    Σs_with_η = copy(S.Σs)

    # Old version: Add -im * η to all elements
    # Σs_with_η .+= -im * η

    # New version: make |Im Σ| at least η
    inds = @. -imag(Σs_with_η) < η
    Σs_with_η[inds] .= real.(Σs_with_η[inds]) .- im * η

    linear_interpolation((S.ωs, S.ks), Σs_with_η; extrapolation_bc = Flat())
end

function get_Σ_itp_dense(S :: PhononDefectSolver, η = S.η)
    Σs_with_η = copy(S.Σs_dense)
    if η != S.η
        # S.Σs_dense is set with S.η, so we need to recompute it.
        Σ_itp_coarse = get_Σ_itp(S, η)
        @views for (ik, k) in enumerate(S.ks_dense)
            @. Σs_with_η[:, ik] = Σ_itp_coarse(S.ωs_dense, k)
        end
    end

    # New version: make |Im Σ| at least η
    inds = @. -imag(Σs_with_η) < η
    Σs_with_η[inds] .= real.(Σs_with_η[inds]) .- im * η

    linear_interpolation((S.ωs_dense, S.ks_dense), Σs_with_η; extrapolation_bc = Flat())
end

"""
    compute_spectral_function!(S :: AbstractSolver; spectral_smearing = S.η)
"""
function compute_phonon_spectral_function!(S :: PhononDefectSolver; spectral_smearing = S.η)
    # Interpolate the self-energy from (ωs, ks) to (ωs_dense, ks_dense)
    Σ_itp = get_Σ_itp(S, spectral_smearing)

    # Compute spectral function
    @views for (ik, k) in enumerate(S.ks_dense)
        #εk = get_εk(k, S.model)
        ωq_sq = get_ωq_sq(k, S.model)
        @. S.Σs_dense[:, ik] = Σ_itp(S.ωs_dense, k)
        @. S.D0[:, ik] = 1 / ((S.ωs_dense .+ 1im * S.η  ).^2 .- ωq_sq )
    #    @. S.As_dense[:, ik] = -imag(1 / (((S.ωs_dense[ik] .+ 1im * S.η ).^2 .- ωq_sq[ik] ) - Σ_itp(S.ωs_dense, k))) / π
       # @. S.As_dense[:, ik] = -imag(1 / ((S.ωs_dense .+ 1im * S.η ).^2 .- ωq_sq - Σ_itp(S.ωs_dense, k) ))/ π
        @. S.As_dense[:, ik] = -imag(S.D0[:, ik] .+  abs.(S.D0[:, ik]).^2  .* S.Σs_dense[:, ik] ) / π

    end

    return S
end


function compute_occupation!(S :: PhononDefectSolver)
    dω = S.ωs_dense[2] - S.ωs_dense[1]
    fermi = @. 1 / (exp((S.ωs_dense - S.model.μ) / S.model.T) + 1)

    @views for (ik, k) in enumerate(S.ks_dense)
        # Compute integral of the spectral function (should be 1 when exact)
        S.spectral_sum[ik] = sum(S.As_dense[:, ik]) * dω

        # Compute occupation
        S.spectral_occ[ik] = sum(S.As_dense[:, ik] .* fermi) * dω
    end

    # Compute total occupation
    if isfinite(S.model.μ)
        # Linearly interpolate the occupations
        occ_itp = linear_interpolation(S.ks_dense, S.spectral_occ)

        # Integrate over the Brillouin zone
        if S.model isa HolsteinLatticeModel || S.model isa PeierlsLatticeModel
            n = quadgk(occ_itp, -π, π)[1] / 2π
        elseif S.model isa FrohlichModel
            n = quadgk(k -> 1 / 2π^2 * k^2 * occ_itp(k), extrema(S.ks_dense)...)[1]
        else
            error("Model not supported")
        end
    else
        n = zero(S.model.μ)
    end

    return n
end


function compute_occupation(S :: PhononDefectSolver, μ; ω_cutoff = -Inf)
    # Same as compute_occupation!, but do not update S.spectral_occ

    dω = step(S.ωs_dense)
    fermi = @. occ_fermion(S.ωs_dense - μ, S.model.T)
    fermi[S.ωs_dense .< ω_cutoff] .= 0

    spectral_occ = (fermi' * S.As_dense)' .* dω

    # Integrate over the Brillouin zone
    if S.model isa HolsteinLatticeModel || S.model isa PeierlsLatticeModel
        # Linearly interpolate the occupations
        occ_itp = linear_interpolation(S.ks_dense, spectral_occ)
        n = quadgk(occ_itp, -π, π)[1] / 2π
    elseif S.model isa FrohlichModel
        # Constant interpolation
        # Integration weights for \int d^3k / (2π)^3 f(k) = \int_a^b dk k^2 f(k) / (2π^2)
        ks_midpoints = vcat(0, (S.ks_dense[1:end-1] .+ S.ks_dense[2:end]) ./ 2, S.ks_dense[end])
        weights_occ = diff(ks_midpoints.^3) ./ 3 ./ 2π^2
        n = sum(weights_occ .* spectral_occ)
    else
        error("Model not supported")
    end

    return n
end

function compute_occupation_MaxwellBoltzmann(S :: PhononDefectSolver)
    # Same as compute_occupation!, but do not update S.spectral_occ

    dω = S.ωs_dense[2] - S.ωs_dense[1]
    fermi = @. exp(-S.ωs_dense / S.model.T)

    spectral_occ = zero(S.spectral_occ)

    @views for (ik, k) in enumerate(S.ks_dense)
        # Compute occupation
        spectral_occ[ik] = sum(S.As_dense[:, ik] .* fermi) * dω
    end

    # Linearly interpolate the occupations
    occ_itp = linear_interpolation(S.ks_dense, spectral_occ)

    # Integrate over the Brillouin zone
    if S.model isa HolsteinLatticeModel || S.model isa PeierlsLatticeModel
        n = quadgk(occ_itp, -π, π)[1] / 2π
    elseif S.model isa FrohlichModel
        n = quadgk(k -> 1 / 2π^2 * k^2 * occ_itp(k), extrema(S.ks_dense)...)[1]
    else
        error("Model not supported")
    end

    return n
end

function plot_phonon_spectral_function!(ax, S :: PhononDefectSolver;
    bare_band = true,
    chemical_potential = true,
    yscale_fac = 1.0,
    kwargs_plot...
    )

    ks = S.ks_dense
    ωs = S.ωs_dense
    dk = ks[2] - ks[1]
    dω = ωs[2] - ωs[1]
    extent = [ks[1] - dk/2, ks[end] + dk/2, (ωs[1] - dω/2) * yscale_fac, (ωs[end] + dω/2) * yscale_fac]

    img = ax.imshow(S.As_dense; origin="lower", extent, aspect="auto", kwargs_plot...)

    if bare_band
        ax.plot(S.ks_dense, sqrt.(get_ωq_sq.(S.ks_dense, S.model)) .* yscale_fac, c="grey", ls="--", lw=1)
    end
    if chemical_potential
        ax.axhline(S.model.μ * yscale_fac, c="r", ls="--", lw=1)
    end

    return img
end

function plot_self_energy!(ax, S :: PhononDefectSolver, term = :real;
    bare_band = true,
    chemical_potential = true,
    yscale_fac = 1.0,
    kwargs_plot...
    )

    Σ_itp = get_Σ_itp(S, 0)

    ks = S.ks_dense
    ωs = S.ωs_dense
    dk = ks[2] - ks[1]
    dω = ωs[2] - ωs[1]
    extent = [ks[1] - dk/2, ks[end] + dk/2, (ωs[1] - dω/2) * yscale_fac, (ωs[end] + dω/2) * yscale_fac]
    Σs = [Σ_itp(ω, k) for ω in ωs, k in ks]

    func = term == :real ? real : imag

    img = ax.imshow(func.(Σs); origin="lower", extent, aspect="auto", kwargs_plot...)

    if bare_band
        ax.plot(S.ks_dense, get_εk.(S.ks_dense, S.model) .* yscale_fac, c="grey", ls="--", lw=1)
    end
    if chemical_potential
        ax.axhline(S.model.μ * yscale_fac, c="r", ls="--", lw=1)
    end

    return img
end