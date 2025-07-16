# Single atomic chain for the phonon dispersion and spectra
using Statistics

struct SingleAtomicChainModel
    M :: Float64
    C :: Float64
    a :: Float64
    m_imp :: Float64
    conc :: Float64
end

Base.Broadcast.broadcastable(model::SingleAtomicChainModel) = Ref(model)
get_ωq_sq(q::Vector , model::SingleAtomicChainModel) = (4 * model.C / model.M) .* sin.(q .* model.a ./ 2).^2
get_ωq_sq(q::Float64 , model::SingleAtomicChainModel) = (4 * model.C / model.M) * sin.(q * model.a / 2)^2
function compute_phonon_self_energy!(S :: PhononDefectSolver{SingleAtomicChainModel})
    (; M, C, a, m_imp, conc) = S.model
    delta_m = m_imp - M
    println(delta_m)
   

    Σ_itp = get_Σ_itp_dense(S, S.η)
    println(size(Σ_itp))

    
    ωq_sq = get_ωq_sq(S.ks, S.model)
    V_mat_R = -delta_m .* S.ωs.^2
    D0_qω = 1.0 ./ ((S.ωs' .+ 1im * S.η ).^2 .- ωq_sq )
    D0_onsite = vec(mean(D0_qω, dims=1))
    T_mat_R = V_mat_R ./ (1.0 .- D0_onsite .* V_mat_R)

    Σs_all = tmap(S.ks) do k
        T_mat_k = zeros(ComplexF64,length(S.ωs))
        for (iω, ω) in enumerate(S.ωs)
            T_mat_k[iω] = T_mat_R[iω]
        end
        conc .* T_mat_k
    end

    @views for ik in eachindex(S.ks)
        S.Σs[:, ik] .= Σs_all[ik] .+= S.Σs_rest[ik]
    end

    return S
end
