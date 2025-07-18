
using .DSpectral: occ_fermion, occ_boson, get_eph_g, get_εk, kramers_kronig

function compute_electron_dressed_phonon_self_energy!(
    S_e::ElectronPhononSolver{PeierlsLatticeModel}, 
    S_ph::PhononDefectSolver{SingleAtomicChainModel}
    )

    e_model = S_e.model
    μ = e_model.μ
    T = e_model.T
    η = S_e.η
    k_grid = S_e.ks_dense
    ε_grid = S_e.ωs_dense
    Σ = S_e.Σs_dense
    VBZ = 1

    q_grid = S_ph.ks_dense
    AqΩ = S_ph.As_dense
    dq = q_grid[2] - q_grid[1]


    Aph_interp = linear_interpolation((S_ph.ωs_dense, S_ph.ks_dense), AqΩ; extrapolation_bc=Flat())


    Σs_imag_all = tmap(k_grid) do k
        Σs_imag_k = zeros(length(ε_grid))
        for (iε, ε) in enumerate(ε_grid)
            sum_val = 0.0
            for q in q_grid
                kq = k + q
                ε_kq = get_εk(kq, e_model)
                ω = ε - ε_kq
                if ω < 0
                    Aph = 0.0
                else
                    Aph = Aph_interp(ω, q)
                end
                f = occ_fermion(ε_kq - μ, T)
                n = occ_boson(ω, T)
                gq = get_eph_g(k, q, e_model)
                g2 = abs2(gq)
                sum_val += g2 * Aph * (1 - f + n) * dq / VBZ
            end
            Σs_imag_k[iε] = -π .* sum_val  
        end
        Σs_imag_k
    end

    # Kramers-Kronig for real part
    @views for (ik, k) in enumerate(k_grid)
        Σs_imag_k = Σs_imag_all[ik]
        Σs_real_k = kramers_kronig(real.(ε_grid), Σs_imag_k)
        @. Σ[:, ik] = Σs_real_k + im * Σs_imag_k
    end

    S_e.Σs = Σ

    return S_e
end