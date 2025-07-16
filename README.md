# DSpectral.jl

A Julia package for spectral analysis in condensed matter physics, focusing on electron-phonon interactions, self-energy calculations, and spectral function computations.

## Features

- **Electron-phonon interactions**: Peierls substitution models and self-consistent calculations
- **Phonon defect calculations**: Single atomic chain models with impurities 
- **Spectral functions**: Compute and visualize A(ω,k) for electrons and phonons
- **k-point grids**: Efficient polynomial spacing grids for 1D, 2D, and 3D integration
- **Kramers-Kronig transformations**: Analytical and numerical implementations
- **Self-energy calculations**: Interpolation and broadening tools
- **Parallel computing**: Multi-threaded calculations using OhMyThreads.jl

## Quick Start

```julia
# Include the package
include("src/DSpectral.jl")
using .DSpectral

# Create a Peierls model
model = PeierlsLatticeModel(
    g = 0.5,    # coupling strength
    ω₀ = 1.0,   # phonon frequency
    t = 1.0,    # hopping parameter
    μ = -2.0,   # chemical potential
    T = 0.1     # temperature
)

# Setup grids
ωs = range(-8.0, 8.0, step=0.02)
ks = range(-π, π, length=513)
qpts = polynomial_grid_1d(π, 512)

# Create and solve
solver = ElectronPhononSolver(model, ωs, ks, qpts; η=0.01, occupation=NaN)
compute_self_energy!(solver)
compute_spectral_function!(solver)

# Visualize (requires PyPlot)
using PyPlot
fig, ax = subplots()
plot_spectral_function!(ax, solver; cmap="viridis")
```

## Documentation

- **[Complete API Documentation](docs/DSpectral_API_Documentation.md)** - Comprehensive guide with examples and usage instructions
- **[Quick Reference](docs/API_Quick_Reference.md)** - Summary of key functions and usage patterns
- **[Examples](examples/)** - Working example scripts demonstrating different features

## Examples

The `examples/` directory contains ready-to-run scripts:

- `peierls_electron.jl` - Electron spectral functions with Peierls substitution
- `single_atomic_chain.jl` - Phonon spectral functions with defects

## Installation

This package is currently in development. Clone the repository and include the source:

```julia
include("path/to/DSpectral.jl/src/DSpectral.jl")
using .DSpectral
```

## Dependencies

Core dependencies:
- `OhMyThreads.jl` - Parallel processing
- `Interpolations.jl` - Grid interpolation
- `Base.Threads` - Threading support

Optional for examples:
- `PyPlot.jl` - Plotting and visualization

## Physics Background

DSpectral.jl implements computational methods for:

1. **Electron-phonon systems**: Calculate how phonon interactions modify electronic properties
2. **Spectral functions**: Compute A(ω,k) = -Im[G(ω,k)]/π to understand quasiparticle behavior
3. **Self-energy**: Calculate Σ(ω,k) corrections due to many-body interactions
4. **Green's functions**: Solve Dyson equations for interacting systems

## Contributing

This is a research package under active development. Contributions, suggestions, and feedback are welcome!

## License

[Add license information]
 
