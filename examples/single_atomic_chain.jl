include("../src/DSpectral.jl")
using .DSpectral
using PyPlot
#using Interpolations: linear_interpolation, Flat

function run_atomic_chain(;
        T, ω₀, λ,
        t = 1.0,
        occupation = 1e-4,
        μ_init = -2.0,
        nq = 128,
        ωs = range(-8., 8., step=0.1),
        η = 0.0,
        plot_spectral = false,
        verbose = false,
        maxiter = 50,
        run_response = false,
        self_consistent = true,
    )
    if !self_consistent && η == 0.
        throw(ArgumentError("η must be provided and nonzero for G0D0 (self_consistent = false)"))
    end

    g = sqrt(λ * ω₀ * t / 2)
    M = 1
    C = 1
    a = 1
    m_imp = 0.5
    conc = 1


    μ = μ_init
    # μ = -Inf
    # occupation = NaN

    ks = range(-π, π, length=nq+1)
    ωs_dense = ωs
    ks_dense = ks

    qpts = DSpectral.Kpoints(Vector(range(-π, π, length=nq+1)[1:end-1]), fill(1 / nq, nq))
    S = PhononDefectSolver(DSpectral.SingleAtomicChainModel(M,C,a,m_imp,conc), ωs, ks, qpts; ks_dense, ωs_dense, η, occupation)
  #  compute_phonon_self_energy!(S)
    compute_phonon_spectral_function!(S)

    return(; S)
end

μ = -Inf
nq = 1024
ωs = range(0., 4., length=1024)

T = 0.0
ω₀ = 1.0
λ = 0.5

η = 1e-2

println("=== ω₀ = $ω₀, λ = $λ, T = $T ===")

S = run_atomic_chain(; ω₀, T, λ, μ_init = μ, nq=nq, ωs=ωs,η, occupation = NaN, self_consistent = false).S


fig, plotaxes = subplots(2,1,figsize=(4, 6))
ax = plotaxes[1]
vmax = 1e2
img = plot_phonon_spectral_function!(bare_band = true,chemical_potential = false,ax, S; cmap="viridis", norm = PyPlot.matplotlib.colors.LogNorm(; vmin=1e-1, vmax))
#img = plot_phonon_spectral_function!(bare_band = false,chemical_potential = false,ax, S; cmap="viridis", vmin = 0, vmax = vmax / 1)
# ax.set_ylim([-5, 4])
#img = plot_phonon_spectral_function!(bare_band = true,chemical_potential = false,ax, S; cmap="viridis")
plotaxes[1].set_ylim(extrema(S.ωs))
# plotaxes[1].set_ylim([-5, -1])
colorbar(img; ax)



dos = dropdims(sum(S.As_dense[:, 1:end-1], dims = 2), dims = 2) .* step(S.ks_dense)
#fermi = @. occ_fermion(S.ωs_dense - S.model.μ, S.model.T)
plotaxes[2].plot(S.ωs_dense, dos)
#plotaxes[2].plot(S.ωs_dense, dos .* fermi)
#plotaxes[2].plot(S.ωs_dense, cumsum(dos .* fermi))
plotaxes[2].set_yscale("log")

show()



    