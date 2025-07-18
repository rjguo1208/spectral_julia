include("../src/DSpectral.jl")
using .DSpectral
using PyPlot

function run_peierls_eph_dressed(;
        T, ω₀, λ,
        t = 1.0,
        occupation = 1e-4,
        μ_init = -0.0,
        nq = 1024,
        nω = 1024,
        η = 0.01
    )

  
    M = 1
    C = 1
    a = 1
    m_imp = 1
    conc = 0.25

    μ = μ_init
    ks = range(-π, π, length=nq+1)
    ωs = range(-8., 8., length=nω+1)
    ωs_dense = ωs
    ks_dense = ks

    qpts = DSpectral.Kpoints(Vector(range(-π, π, length=nq+1)[1:end-1]), fill(1 / nq, nq))
    S_ph = PhononDefectSolver(DSpectral.SingleAtomicChainModel(M,C,a,m_imp,conc), ωs, ks, qpts; ks_dense, ωs_dense, η, occupation)
    compute_phonon_self_energy!(S_ph)
    compute_phonon_spectral_function!(S_ph)

    # Electron
    g = sqrt(λ * ω₀ * t / 2)
    S_e = ElectronPhononSolver(DSpectral.PeierlsLatticeModel(g, ω₀, t, μ, T), ωs, ks, qpts; ks_dense, ωs_dense, η, occupation)
    compute_electron_dressed_phonon_self_energy!(S_e, S_ph)
    compute_spectral_function!(S_e)

    return(; S_e, S_ph)
end

μ = -2.0
nq = 512
nω = 512

T = 0.0
ω₀ = 1.0
λ = 0.5

η = 1e-2

println("=== ω₀ = $ω₀, λ = $λ, T = $T ===")

S_e, S_ph = run_peierls_eph_dressed(; ω₀, T, λ, μ_init = μ, nq=nq, nω=nω, η, occupation = NaN)
println(S_e.Σs_dense[1, 1])


fig, plotaxes = subplots(figsize=(4, 3))
ax = plotaxes
vmax = 1e2
img = plot_phonon_spectral_function!(bare_band = true,chemical_potential = false,ax, S_ph; cmap="viridis", norm = PyPlot.matplotlib.colors.LogNorm(; vmin=1e-1, vmax))
#img = plot_phonon_spectral_function!(bare_band = false,chemical_potential = false,ax, S; cmap="viridis", vmin = 0, vmax = vmax / 1)
# ax.set_ylim([-5, 4])
#img = plot_phonon_spectral_function!(bare_band = true,chemical_potential = false,ax, S; cmap="viridis")
plotaxes.set_ylim(extrema(S_ph.ωs))
# plotaxes[1].set_ylim([-5, -1])
colorbar(img; ax)
show()

fig, plotaxes = subplots(figsize=(4, 3))
ax = plotaxes
vmax = 10
img = plot_spectral_function!(bare_band = true,chemical_potential = true,ax, S_e; cmap="viridis", norm = PyPlot.matplotlib.colors.LogNorm(; vmin=vmax/1e3, vmax))
#img = plot_phonon_spectral_function!(bare_band = false,chemical_potential = false,ax, S; cmap="viridis", vmin = 0, vmax = vmax / 1)
# ax.set_ylim([-5, 4])
#img = plot_phonon_spectral_function!(bare_band = true,chemical_potential = false,ax, S; cmap="viridis")
plotaxes.set_ylim(extrema(S_e.ωs))
# plotaxes[1].set_ylim([-5, -1])
colorbar(img; ax)
show()




