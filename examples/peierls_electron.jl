include("../src/DSpectral.jl")
using .DSpectral
using PyPlot

function run_peierls(;
        T, ω₀, λ,
        t = 1.0,
        occupation = 1e-4,
        μ_init = -2.0,
        nq = 128,
        ωs = range(-8., 8., step=0.01),
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

    μ = μ_init
    # μ = -Inf
    # occupation = NaN

    ks = range(-π, π, length=nq+1)
    ωs_dense = ωs
    ks_dense = ks

    qpts = DSpectral.Kpoints(Vector(range(-π, π, length=nq+1)[1:end-1]), fill(1 / nq, nq))
    S = ElectronPhononSolver(DSpectral.PeierlsLatticeModel(g, ω₀, t, μ, T), ωs, ks, qpts; ks_dense, ωs_dense, η, occupation)
    
    if self_consistent  # scGD0
        ;
       # EPSpectral.ep_solve!(S; tol=1e-5, maxiter, update_μ = !isnan(occupation), verbose);
    else  # G0D0
        DSpectral.compute_self_energy!(S)
    end
    compute_spectral_function!(S)
    return (; S)
end


μ = -Inf
nq = 512
ωs = range(-8., 8., step=0.02)

T = 0.0
ω₀ = 1.0
λ = 0.5

η = 1e-2

println("=== ω₀ = $ω₀, λ = $λ, T = $T ===")

S = run_peierls(; ω₀, T, λ, μ_init = μ, nq, η, occupation = NaN, self_consistent = false).S

fig, plotaxes = subplots(2, 1, figsize=(5, 8))
ax = plotaxes[1]
vmax = maximum(S.As_dense)
img = plot_spectral_function!(ax, S; cmap="viridis", norm = PyPlot.matplotlib.colors.LogNorm(; vmin=vmax/1e3, vmax))
# img = plot_spectral_function!(ax, S; cmap="viridis", vmin = 0, vmax = vmax / 1)
# ax.set_ylim([-5, 4])
plotaxes[1].set_ylim(extrema(S.ωs))
# plotaxes[1].set_ylim([-5, -1])
colorbar(img; ax)

dos = dropdims(sum(S.As_dense[:, 1:end-1], dims = 2), dims = 2) .* step(S.ks_dense)
fermi = @. occ_fermion(S.ωs_dense - S.model.μ, S.model.T)
plotaxes[2].plot(S.ωs_dense, dos)
plotaxes[2].plot(S.ωs_dense, dos .* fermi)
plotaxes[2].plot(S.ωs_dense, cumsum(dos .* fermi))
plotaxes[2].set_yscale("log")
show()
