# DSpectral.jl - Quick API Reference

## Core Functions (Exported)

### Grid Generation
```julia
polynomial_grid_1d(kmax, n, order=1)      # 1D polynomial spacing grid
polynomial_grid_2d(kmax, n, order=1)      # 2D grid inside sphere
polynomial_grid_3d(kmax, n, order=1)      # 3D grid inside sphere  
azimuthal_grid_3d(kmax, nk, nθ, order=1)  # 3D azimuthally symmetric grid
```

### Models
```julia
PeierlsLatticeModel(g, ω₀, t, μ, T)           # Electron-phonon model
SingleAtomicChainModel(M, C, a, m_imp, conc)  # Phonon chain with defects
```

### Solvers
```julia
ElectronPhononSolver(model, ωs, ks, qpts; kwargs...)  # Electron solver
PhononDefectSolver(model, ωs, ks, qpts; kwargs...)    # Phonon solver
```

### Self-Energy & Spectral Functions
```julia
get_Σ_itp(solver, η)                          # Self-energy interpolation
get_Σ_itp_dense(solver, η)                    # Dense grid interpolation
compute_spectral_function!(solver)             # Electron spectral function
compute_phonon_spectral_function!(solver)     # Phonon spectral function
compute_phonon_self_energy!(solver)           # Phonon self-energy
```

### Plotting
```julia
plot_spectral_function!(ax, solver; kwargs...)        # Plot electron A(ω,k)
plot_phonon_spectral_function!(ax, solver; kwargs...) # Plot phonon A(ω,k)
```

### Utility Functions
```julia
occ_fermion(e, T)                # Fermi-Dirac distribution
occ_boson(e, T)                  # Bose-Einstein distribution  
delta_gaussian(x, η)             # Gaussian δ-function
delta_lorentzian(x, η)           # Lorentzian δ-function
kramers_kronig(ωs, zs; tail)     # Kramers-Kronig transform
```

## Quick Usage Patterns

### Basic Electron-Phonon Calculation
```julia
# Setup
model = PeierlsLatticeModel(g, ω₀, t, μ, T)
solver = ElectronPhononSolver(model, ωs, ks, qpts; η=η, occupation=occ)

# Solve
compute_self_energy!(solver)  
compute_spectral_function!(solver)

# Plot
fig, ax = subplots()
plot_spectral_function!(ax, solver)
```

### Basic Phonon Defect Calculation  
```julia
# Setup
model = SingleAtomicChainModel(M, C, a, m_imp, conc)
solver = PhononDefectSolver(model, ωs, ks, qpts; η=η, occupation=NaN)

# Solve
compute_phonon_self_energy!(solver)
compute_phonon_spectral_function!(solver)

# Plot
fig, ax = subplots()
plot_phonon_spectral_function!(ax, solver)
```

## Key Data Access

### From ElectronPhononSolver
```julia
solver.As_dense        # Electron spectral function A(ω,k)
solver.Σs_dense        # Self-energy Σ(ω,k) on dense grid
solver.ωs_dense        # Dense frequency grid
solver.ks_dense        # Dense k-point grid
solver.model           # Physics model
solver.η               # Broadening parameter
```

### From PhononDefectSolver
```julia
solver.As_dense        # Phonon spectral function A(ω,k)
solver.Σs_dense        # Phonon self-energy on dense grid
solver.D0              # Bare phonon Green's function
solver.ωs_dense        # Dense frequency grid
solver.ks_dense        # Dense k-point grid
```

## Common Workflows

### 1. Parameter Scan
```julia
λ_values = [0.1, 0.2, 0.5, 1.0]
results = []

for λ in λ_values
    g = sqrt(λ * ω₀ * t / 2)
    model = PeierlsLatticeModel(g, ω₀, t, μ, T)
    solver = ElectronPhononSolver(model, ωs, ks, qpts; η=η, occupation=occ)
    
    compute_self_energy!(solver)
    compute_spectral_function!(solver)
    
    push!(results, (λ=λ, solver=solver))
end
```

### 2. Self-Consistent Loop
```julia
for iter in 1:maxiter
    Σ_old = copy(solver.Σs)
    compute_self_energy!(solver)
    
    error = maximum(abs.(solver.Σs - Σ_old))
    if error < tolerance
        break
    end
end
```

### 3. Extract Physical Quantities
```julia
# Density of states
dos = dropdims(sum(solver.As_dense[:, 1:end-1], dims=2), dims=2) * step(solver.ks_dense)

# Spectral weight per k-point  
spectral_weight = sum(solver.As_dense, dims=1) * step(solver.ωs_dense)

# Interpolate self-energy at specific point
Σ_itp = get_Σ_itp_dense(solver)
Σ_value = Σ_itp(ω_test, k_test)
```

## Typical Parameter Ranges

### Peierls Model
- `g`: 0.1 - 2.0 (coupling strength)
- `ω₀`: 0.5 - 2.0 (phonon frequency)  
- `t`: 1.0 (hopping, often used as energy unit)
- `μ`: -3.0 to 3.0 (chemical potential)
- `T`: 0.01 - 1.0 (temperature)
- `η`: 1e-3 - 1e-1 (broadening)

### Grid Sizes
- `nω`: 100 - 2000 (frequency points)
- `nk`: 50 - 1000 (k-points) 
- `nq`: 100 - 1000 (q-points for integration)
- Dense grids: 3-10x denser than coarse grids

## Performance Notes

- Use `polynomial_grid_*` functions for efficient integration
- Dense grids improve spectral function smoothness
- Increase `η` if encountering numerical instabilities  
- Package uses threading via `OhMyThreads` - set `JULIA_NUM_THREADS` appropriately
- Memory usage scales as `nω × nk` for spectral functions

## See Also

- [Complete API Documentation](DSpectral_API_Documentation.md) - Full documentation with examples
- [Examples](../examples/) - Working example scripts
- [Source Code](../src/) - Implementation details