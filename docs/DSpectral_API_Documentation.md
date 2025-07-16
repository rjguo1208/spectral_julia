# DSpectral.jl - Comprehensive API Documentation

DSpectral.jl is a Julia package for spectral analysis in condensed matter physics, focusing on electron-phonon interactions, self-energy calculations, and spectral function computations.

## Table of Contents

1. [Installation and Setup](#installation-and-setup)
2. [Core Data Types](#core-data-types)
3. [Grid Generation Functions](#grid-generation-functions)
4. [Physics Models](#physics-models)
5. [Solver Types](#solver-types)
6. [Self-Energy and Spectral Functions](#self-energy-and-spectral-functions)
7. [Utility Functions](#utility-functions)
8. [Kramers-Kronig Transformations](#kramers-kronig-transformations)
9. [Complete Examples](#complete-examples)
10. [Advanced Features](#advanced-features)

## Installation and Setup

```julia
# Include the package (assuming local development)
include("src/DSpectral.jl")
using .DSpectral

# For plotting examples (optional)
using PyPlot
```

## Core Data Types

### Kpoints{T}

A structure for representing k-point grids with associated integration weights.

**Fields:**
- `vectors :: Vector{T}`: k-point vectors
- `weights :: Vector{Float64}`: Integration weights

**Interface:**
```julia
# Indexing and iteration
Base.length(k::Kpoints) = length(k.vectors)
Base.getindex(k::Kpoints, i::Int) = (k.vectors[i], k.weights[i])
Base.iterate(k::Kpoints, i=1) = i > length(k) ? nothing : (k[i], i+1)

# Get dimension
get_dimension(k::Kpoints) = length(k.vectors[1])
```

**Example:**
```julia
# Create a simple k-point grid
vectors = [0.1, 0.2, 0.3, 0.4]
weights = [0.25, 0.25, 0.25, 0.25]
kpts = DSpectral.Kpoints(vectors, weights)

# Iterate through k-points
for (k, w) in kpts
    println("k = $k, weight = $w")
end
```

## Grid Generation Functions

### polynomial_grid_1d(kmax, n, order=1)

Generate a 1D grid with polynomial spacing in [-kmax, kmax].

**Parameters:**
- `kmax`: Maximum k-value
- `n :: Int`: Number of grid points
- `order :: Int = 1`: Polynomial order for spacing

**Returns:** `Kpoints` object with 1D grid

**Examples:**
```julia
# Linear spacing (order=1)
kpts_1d = polynomial_grid_1d(π, 10, 1)

# Quadratic spacing (order=2) - denser near boundaries
kpts_1d_quad = polynomial_grid_1d(π, 10, 2)

# Access grid points and weights
for (k, w) in kpts_1d
    println("k = $k, weight = $w")
end
```

### polynomial_grid_2d(kmax, n, order=1)

Generate a 2D grid with polynomial spacing inside a sphere of radius kmax.

**Parameters:**
- `kmax`: Maximum radius
- `n :: Int`: Number of grid points per dimension
- `order :: Int = 1`: Polynomial order

**Returns:** `Kpoints` object with 2D grid (filtered for |k| ≤ kmax)

**Example:**
```julia
# 2D grid for Brillouin zone integration
kpts_2d = polynomial_grid_2d(π, 20, 1)
println("Total k-points in 2D grid: $(length(kpts_2d))")
```

### polynomial_grid_3d(kmax, n, order=1)

Generate a 3D grid with polynomial spacing inside a sphere of radius kmax.

**Parameters:**
- `kmax`: Maximum radius  
- `n :: Int`: Number of grid points per dimension
- `order :: Int = 1`: Polynomial order

**Returns:** `Kpoints` object with 3D grid (filtered for |k| ≤ kmax)

**Example:**
```julia
# 3D grid for bulk calculations
kpts_3d = polynomial_grid_3d(2π, 15, 1)
println("Total k-points in 3D grid: $(length(kpts_3d))")
```

### azimuthal_grid_3d(kmax, nk, nθ, order=1; kmin=0.0)

Generate a 3D grid with azimuthal symmetry along the x-axis.

**Parameters:**
- `kmax`: Maximum k-value
- `nk :: Int`: Number of radial points
- `nθ :: Int`: Number of angular points
- `order :: Int = 1`: Polynomial order for radial grid
- `kmin = 0.0`: Minimum k-value

**Returns:** `Kpoints` object with azimuthally symmetric grid

**Example:**
```julia
# Grid for cylindrically symmetric problems
kpts_cyl = azimuthal_grid_3d(π, 50, 20, 1)
```

## Physics Models

### PeierlsLatticeModel

Model for electron-phonon interaction using Peierls substitution.

**Fields:**
- `g :: Float64`: Electron-phonon coupling strength
- `ω₀ :: Float64`: Phonon frequency
- `t :: Float64`: Hopping parameter
- `μ :: Float64`: Chemical potential
- `T :: Float64`: Temperature

**Key Functions:**
```julia
# Dispersion relation
get_εk(k, model::PeierlsLatticeModel) = -2 * model.t * cos(k)

# Velocity
get_vk(k, model::PeierlsLatticeModel) = 2 * model.t * sin(k)

# Electron-phonon coupling
get_eph_g(k, q, model::PeierlsLatticeModel) = -2 * im * model.g * (sin(k + q) - sin(k))
```

**Example:**
```julia
# Create Peierls model
model = PeierlsLatticeModel(
    g = 0.5,    # coupling strength
    ω₀ = 1.0,   # phonon frequency 
    t = 1.0,    # hopping
    μ = -2.0,   # chemical potential
    T = 0.1     # temperature
)

# Calculate dispersion at k=π/4
k = π/4
ε = get_εk(k, model)
v = get_vk(k, model)
println("ε(k=$k) = $ε, v(k=$k) = $v")
```

### SingleAtomicChainModel

Model for phonon dispersion in a single atomic chain with impurities.

**Fields:**
- `M :: Float64`: Atomic mass
- `C :: Float64`: Spring constant
- `a :: Float64`: Lattice spacing
- `m_imp :: Float64`: Impurity mass
- `conc :: Float64`: Impurity concentration

**Key Functions:**
```julia
# Phonon dispersion (squared frequency)
get_ωq_sq(q, model::SingleAtomicChainModel) = (4 * model.C / model.M) * sin(q * model.a / 2)^2
```

**Example:**
```julia
# Create atomic chain model
chain_model = SingleAtomicChainModel(
    M = 1.0,      # host mass
    C = 1.0,      # spring constant
    a = 1.0,      # lattice spacing
    m_imp = 0.5,  # impurity mass
    conc = 0.1    # 10% impurity concentration
)

# Calculate phonon frequency
q = π/2
ω²= get_ωq_sq(q, chain_model)
ω = sqrt(ω²)
println("ω(q=$q) = $ω")
```

## Solver Types

### ElectronPhononSolver{MT, VT}

Main solver for electron-phonon problems with self-consistent calculations.

**Constructor:**
```julia
ElectronPhononSolver(model, ωs, ks, qpts; 
                    occupation, ωs_dense=nothing, ks_dense=nothing, 
                    η, Σs_rest=nothing)
```

**Parameters:**
- `model`: Physics model (e.g., `PeierlsLatticeModel`)
- `ωs`: Frequency grid (can be non-uniform)
- `ks`: k-point grid (can be non-uniform)  
- `qpts`: q-point grid for integrations
- `occupation`: Target electron occupation
- `ωs_dense`, `ks_dense`: Dense grids for interpolation
- `η`: Broadening parameter
- `Σs_rest`: Additional self-energy contributions

**Example:**
```julia
# Setup grids
ωs = range(-8.0, 8.0, step=0.01)
ks = range(-π, π, length=129)
nq = 128
qpts = DSpectral.Kpoints(
    Vector(range(-π, π, length=nq+1)[1:end-1]), 
    fill(1/nq, nq)
)

# Create solver
model = PeierlsLatticeModel(0.5, 1.0, 1.0, -2.0, 0.1)
solver = ElectronPhononSolver(
    model, ωs, ks, qpts;
    occupation = 1e-4,
    η = 0.01
)

# Compute self-energy and spectral function
compute_self_energy!(solver)
compute_spectral_function!(solver)
```

### PhononDefectSolver{MT, VT}

Solver for phonon problems with defects and impurities.

**Constructor:**
```julia
PhononDefectSolver(model, ωs, ks, qpts; 
                  occupation, ωs_dense=nothing, ks_dense=nothing,
                  η, Σs_rest=nothing)
```

**Example:**
```julia
# Setup for phonon calculation
ωs = range(0.0, 4.0, step=0.005)
ks = range(-π, π, length=1025)
qpts = DSpectral.Kpoints(Vector(ks[1:end-1]), fill(1/1024, 1024))

# Create phonon solver
chain_model = SingleAtomicChainModel(1.0, 1.0, 1.0, 0.5, 0.1)
phonon_solver = PhononDefectSolver(
    chain_model, ωs, ks, qpts;
    occupation = NaN,
    η = 1e-3
)

# Compute phonon self-energy and spectral function
compute_phonon_self_energy!(phonon_solver)
compute_phonon_spectral_function!(phonon_solver)
```

## Self-Energy and Spectral Functions

### get_Σ_itp(S::AbstractSolver, η=S.η)

Create interpolation object for self-energy with specified broadening.

**Parameters:**
- `S`: Solver object
- `η`: Broadening parameter (defaults to solver's η)

**Returns:** Interpolation object for Σ(ω, k)

### get_Σ_itp_dense(S::AbstractSolver, η=S.η) 

Create interpolation object for self-energy on dense grid.

### compute_spectral_function!(S::ElectronPhononSolver; spectral_smearing=S.η)

Compute the electron spectral function A(ω,k) = -Im[G(ω,k)]/π.

**Example:**
```julia
# After creating and solving
compute_spectral_function!(solver; spectral_smearing=0.01)

# Access spectral function
A_ωk = solver.As_dense  # Matrix: A(ω,k)
```

### compute_phonon_spectral_function!(S::PhononDefectSolver; spectral_smearing=S.η)

Compute the phonon spectral function.

### plot_spectral_function!(ax, S::ElectronPhononSolver; kwargs...)

Plot the electron spectral function as a 2D colormap.

**Example:**
```julia
using PyPlot
fig, ax = subplots()
img = plot_spectral_function!(ax, solver; cmap="viridis")
colorbar(img; ax)
ax.set_xlabel("k")
ax.set_ylabel("ω")
```

### plot_phonon_spectral_function!(ax, S::PhononDefectSolver; bare_band=true, chemical_potential=false, kwargs...)

Plot the phonon spectral function.

## Utility Functions

### occ_fermion(e, T)

Fermi-Dirac distribution function.

**Parameters:**
- `e`: Energy (can be relative to chemical potential)
- `T`: Temperature

**Returns:** Occupation probability f(e) = 1/(exp(e/T) + 1)

**Example:**
```julia
# Calculate occupation at different temperatures
energies = range(-3, 3, length=100)
T_cold = 0.1
T_hot = 1.0

occ_cold = occ_fermion.(energies, T_cold)
occ_hot = occ_fermion.(energies, T_hot)
```

### occ_boson(e, T)

Bose-Einstein distribution function.

**Parameters:**
- `e`: Energy
- `T`: Temperature  

**Returns:** Occupation n(e) = 1/(exp(e/T) - 1)

### delta_gaussian(x, η)

Gaussian representation of delta function.

**Parameters:**
- `x`: Argument
- `η`: Broadening parameter

**Returns:** δ(x) ≈ exp(-(x/η)²)/(√π η)

### delta_lorentzian(x, η)

Lorentzian representation of delta function.

**Parameters:**
- `x`: Argument  
- `η`: Broadening parameter

**Returns:** δ(x) ≈ η/(π(x² + η²))

**Example:**
```julia
# Compare delta function representations
x = range(-1, 1, length=1000)
η = 0.1

δ_gauss = delta_gaussian.(x, η)
δ_lorentz = delta_lorentzian.(x, η)

# Should integrate to approximately 1
∫_gauss = trapz(x, δ_gauss)   # ≈ 1
∫_lorentz = trapz(x, δ_lorentz) # ≈ 1
```

## Kramers-Kronig Transformations

### kramers_kronig(ωs, zs; tail=false)

Calculate the Kramers-Kronig transformation to obtain the real part from the imaginary part.

**Parameters:**
- `ωs`: Frequency grid
- `zs`: Imaginary part of retarded function
- `tail`: Whether to include tail corrections

**Returns:** Real part ys such that ys + i*zs is analytic in upper half-plane

**Mathematical Formula:**
y(ω) = (1/π) ∫ dω' z(ω')/(ω' - ω)

**Example:**
```julia
# Example: Kramers-Kronig transform of a Lorentzian
ωs = range(-5, 5, length=1000)
η = 0.1
zs = -η ./ (ωs.^2 .+ η^2)  # Lorentzian imaginary part

# Get real part via KK transform
ys = kramers_kronig(ωs, zs)

# The result should be: ys ≈ ωs ./ (ωs.^2 .+ η^2)
```

## Complete Examples

### Example 1: Peierls Model Electron Spectral Function

```julia
include("src/DSpectral.jl")
using .DSpectral
using PyPlot

function peierls_example()
    # Parameters
    g = 0.5      # coupling
    ω₀ = 1.0     # phonon frequency
    t = 1.0      # hopping
    μ = -2.0     # chemical potential 
    T = 0.1      # temperature
    η = 0.01     # broadening
    
    # Grids
    ωs = range(-8.0, 8.0, step=0.02)
    nq = 512
    ks = range(-π, π, length=nq+1)
    qpts = DSpectral.Kpoints(
        Vector(ks[1:end-1]), 
        fill(1/nq, nq)
    )
    
    # Create and solve
    model = PeierlsLatticeModel(g, ω₀, t, μ, T)
    solver = ElectronPhononSolver(
        model, ωs, ks, qpts;
        ks_dense=ks, ωs_dense=ωs,
        η=η, occupation=NaN
    )
    
    compute_self_energy!(solver)
    compute_spectral_function!(solver)
    
    # Plot results
    fig, (ax1, ax2) = subplots(2, 1, figsize=(6, 8))
    
    # Spectral function
    vmax = maximum(solver.As_dense)
    img = plot_spectral_function!(ax1, solver; 
        cmap="viridis", 
        norm=PyPlot.matplotlib.colors.LogNorm(vmin=vmax/1e3, vmax=vmax)
    )
    colorbar(img; ax=ax1)
    ax1.set_title("Electron Spectral Function A(ω,k)")
    
    # Density of states
    dos = dropdims(sum(solver.As_dense[:, 1:end-1], dims=2), dims=2) * step(solver.ks_dense)
    ax2.plot(solver.ωs_dense, dos)
    ax2.set_xlabel("ω")
    ax2.set_ylabel("DOS")
    ax2.set_yscale("log")
    ax2.set_title("Density of States")
    
    tight_layout()
    show()
    
    return solver
end

# Run the example
solver = peierls_example()
```

### Example 2: Phonon Spectral Function with Defects

```julia
function phonon_defect_example()
    # Parameters
    M = 1.0       # host mass
    C = 1.0       # spring constant
    a = 1.0       # lattice spacing
    m_imp = 0.5   # impurity mass
    conc = 0.1    # concentration
    η = 1e-3      # broadening
    
    # Grids
    ωs = range(0.0, 4.0, step=0.005)
    nq = 1024
    ks = range(-π, π, length=nq+1)
    qpts = DSpectral.Kpoints(
        Vector(ks[1:end-1]),
        fill(1/nq, nq)
    )
    
    # Create and solve
    model = SingleAtomicChainModel(M, C, a, m_imp, conc)
    solver = PhononDefectSolver(
        model, ωs, ks, qpts;
        ks_dense=ks, ωs_dense=ωs,
        η=η, occupation=NaN
    )
    
    compute_phonon_self_energy!(solver)
    compute_phonon_spectral_function!(solver)
    
    # Plot
    fig, ax = subplots(figsize=(6, 4))
    img = plot_phonon_spectral_function!(ax, solver;
        bare_band=true, chemical_potential=false,
        cmap="viridis",
        norm=PyPlot.matplotlib.colors.LogNorm(vmin=1e-1, vmax=1e2)
    )
    colorbar(img; ax=ax)
    ax.set_title("Phonon Spectral Function with Defects")
    show()
    
    return solver
end

# Run the phonon example
phonon_solver = phonon_defect_example()
```

## Advanced Features

### Self-Consistent Calculations

For self-consistent calculations, iterate the solver until convergence:

```julia
function self_consistent_solve!(solver; tol=1e-5, maxiter=50)
    for iter in 1:maxiter
        Σ_old = copy(solver.Σs)
        
        compute_self_energy!(solver)
        compute_spectral_function!(solver)
        
        # Check convergence
        error = maximum(abs.(solver.Σs - Σ_old))
        println("Iteration $iter: error = $error")
        
        if error < tol
            println("Converged!")
            break
        end
    end
    return solver
end
```

### Custom Interpolation and Analysis

```julia
# Extract and analyze results
function analyze_spectral_function(solver)
    # Get interpolation object
    Σ_itp = get_Σ_itp_dense(solver)
    
    # Calculate self-energy at specific points
    ω_test, k_test = 0.0, π/4
    Σ_value = Σ_itp(ω_test, k_test)
    println("Σ(ω=$ω_test, k=$k_test) = $Σ_value")
    
    # Calculate spectral weight
    dω = step(solver.ωs_dense)
    dk = step(solver.ks_dense)
    
    total_weight = sum(solver.As_dense) * dω * dk / (2π)
    println("Total spectral weight: $total_weight")
    
    return Σ_itp
end
```

### Performance Tips

1. **Grid Selection**: Use polynomial grids with appropriate order for your problem
2. **Broadening**: Choose η to balance accuracy and stability
3. **Dense Grids**: Use dense grids for smooth spectral functions
4. **Parallel Processing**: The package uses `OhMyThreads` for parallelization

### Troubleshooting

**Common Issues:**

1. **Convergence Problems**: 
   - Increase broadening parameter η
   - Reduce grid spacing
   - Check model parameters

2. **Memory Issues**:
   - Reduce dense grid size
   - Use smaller frequency/k-point ranges

3. **Numerical Instabilities**:
   - Ensure proper broadening
   - Check for singular points in dispersion

**Debug Information:**
```julia
# Print solver information
println(solver)  # Shows grid sizes and parameters

# Check spectral sum rule
spectral_sum = sum(solver.As_dense, dims=1) * step(solver.ωs_dense)
println("Spectral sum per k-point: ", extrema(spectral_sum))
```

This documentation covers the complete public API of DSpectral.jl with practical examples and usage instructions for condensed matter physics calculations involving electron-phonon interactions and spectral functions.