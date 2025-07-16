
struct Kpoints{T}
    vectors :: Vector{T}
    weights :: Vector{Float64}
end

# Indexing and iteration interface
Base.length(k :: Kpoints) = length(k.vectors)
Base.getindex(k :: Kpoints, i::Int) = (k.vectors[i], k.weights[i])
Base.firstindex(k :: Kpoints) = 1
Base.lastindex(k :: Kpoints) = length(k)
Base.iterate(k :: Kpoints, i = 1) = i > length(k) ? nothing : (k[i], i + 1)

get_dimension(k :: Kpoints) = length(k.vectors[1])



"""
    uniform_grid_1d(kmax, n)
Generate a 1D grid with uniform spacing in [-kmax,kmax].
"""
# TODO

"""
    polynomial_grid_1d(kmax, n, order = 1)
Generate a 1D grid with polynomial spacing in [-kmax, kmax].
For n = 2m - 1: ``k[i+m] = kmax * (i / (m - 1))^order``  (i = 1, ..., m - 1)
For n = 2m    : ``k[i+m] = kmax * ((i - 1/2) / (m - 1/2))^order`` (i = 1, ..., m)
"""
function polynomial_grid_1d(kmax, n :: Int, order :: Int = 1)
    if mod(n, 2) == 1
        m = div(n + 1, 2)
        ks_positive = @. kmax * ((1:(m-1)) ./ (m - 1))^order
        ks = vcat(.-reverse(ks_positive), 0, ks_positive)
    else
        m = div(n, 2)
        ks_positive = @. kmax * (((1:m) - 1/2) ./ (m - 1/2))^order
        ks = vcat(.-reverse(ks_positive), ks_positive)
    end

    weights_ = (ks[3:end] .- ks[1:end-2]) ./ 2
    weights = vcat(weights_[1], weights_, weights_[end])

    Kpoints(ks, weights)
end

"""
    polynomial_grid_2d(kmax, n, order = 1)
Generate a 2D grid with polynomial spacing in each directions inside the sphere with radius
kmax.
"""
function polynomial_grid_2d(kmax, n :: Int, order :: Int = 1)
    kpts_1d = polynomial_grid_1d(kmax, n, order)
    vectors_1d, weights_1d = kpts_1d.vectors, kpts_1d.weights

    vectors = vec(SVector.(collect(Iterators.product(vectors_1d, vectors_1d))))
    weights = vec(prod.(collect(Iterators.product(weights_1d, weights_1d))))

    # Filter out k points with |k| > kmax
    inds = norm.(vectors) .<= kmax
    Kpoints(vectors[inds], weights[inds])
end


"""
    polynomial_grid_3d(kmax, n, order = 1)
Generate a 3D grid with polynomial spacing in each directions inside the sphere with radius
kmax.
"""
function polynomial_grid_3d(kmax, n :: Int, order :: Int = 1)
    kpts_1d = polynomial_grid_1d(kmax, n, order)
    vectors_1d, weights_1d = kpts_1d.vectors, kpts_1d.weights

    vectors = vec(SVector.(collect(Iterators.product(vectors_1d, vectors_1d, vectors_1d))))
    weights = vec(prod.(collect(Iterators.product(weights_1d, weights_1d, weights_1d))))

    # Filter out k points with |k| > kmax
    inds = norm.(vectors) .<= kmax
    Kpoints(vectors[inds], weights[inds])
end


"""
    azimuthal_grid_3d(kmax, n, order)
Generate a 3D grid with azimuthal symmetry along the x axis.
"""
function azimuthal_grid_3d(kmax, nk, nθ, order = 1; kmin = 0.0)
    # Radial grid: k² dk
    if kmin != 0
        order != 1 && throw(ArgumentError("order must be 1 for kmin != 0"))
        ks = range(kmin, kmax, length = nk)
    else
        ks = @. kmax * ((1:nk) ./ nk)^order
    end

    weights_k = zeros(nk)
    for ik in 1:nk
        if ik < nk
            a, b = ks[ik], ks[ik+1]
            weights_k[ik] += (b - a) * (3*a^2 + 2*a*b + b^2) / 12
        end
        if ik > 1
            a, b = ks[ik-1], ks[ik]
            weights_k[ik] += (b - a) * (a^2 + 2*a*b + 3*b^2) / 12
        end
    end

    # Angular grid: sinθ dθ = d(cosθ)
    # Compute weights from the Gauss-Legendre quadrature
    cosθs, weights_θ = gausslegendre(nθ)

    vectors = map(Iterators.product(ks, cosθs)) do (k, cosθ)
        sinθ = sqrt(1 - cosθ^2)
        SVector(k * cosθ, k * sinθ, 0)
    end |> vec

    # Multiply 2π for the azimuthal weight
    weights = vec(prod.(collect(Iterators.product(weights_k, weights_θ)))) .* 2π

    Kpoints(vectors, weights)
end

"""
    spherical_grid_3d(kmax, n, order)
Generate a 3D grid assuming a spherically symmetric integrand.
"""
function spherical_grid_3d(kmax, nk, order = 1)
    # Radial grid: k² dk
    ks = @. kmax * ((1:nk) ./ nk)^order

    weights_k = zeros(nk)
    for ik in 1:nk
        if ik < nk
            a, b = ks[ik], ks[ik+1]
            weights_k[ik] += (b - a) * (3*a^2 + 2*a*b + b^2) / 12
        end
        if ik > 1
            a, b = ks[ik-1], ks[ik]
            weights_k[ik] += (b - a) * (a^2 + 2*a*b + 3*b^2) / 12
        end
    end

    vectors = SVector.(ks, 0, 0)

    # Multiply 4π for the spherical weight
    weights = weights_k .* 4π

    Kpoints(vectors, weights)
end
