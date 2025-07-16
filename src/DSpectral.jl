__precompile__(true)
module DSpectral
    using PrecompileTools

    @recompile_invalidations begin
       # using LinearAlgebra
       # using StaticArrays
        using OhMyThreads
       # using ChunkSplitters
       # using PolyLog
        using Interpolations
       # using QuadGK
       # using Roots
       # using NLsolve
       # using FastGaussQuadrature
        using Base.Threads
    end
    include("utils.jl")
    include("kpoints.jl")
    include("kramers_kronig.jl")
    include("solver.jl")
    include("peierls.jl")
    include("linear_spline_basis.jl")
    include("selfen.jl")
    include("solver_pd.jl")
    include("atomic_chain.jl")
    include("selfen_phonon.jl")


    export
    polynomial_grid_1d, polynomial_grid_2d, polynomial_grid_3d, azimuthal_grid_3d,
    kramers_kronig, occ_fermion,
    PeierlsLatticeModel,SingleAtomicChainModel,
    ElectronPhononSolver,PhononDefectSolver,
    get_Σ_itp, get_Σ_itp_dense,
    compute_spectral_function!,
    plot_spectral_function!,
    plot_phonon_spectral_function!,
    compute_phonon_self_energy!,
    compute_phonon_spectral_function!
end
