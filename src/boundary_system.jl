struct BoundarySystemM
    A::Matrix{CT}
    B::Matrix{CT}
    weights::Vector{Float64}
    m::Int
    nmin::Int
    nmax::Int
end

function _check_axisymmetric_shape(shape)
    hasmethod(TransitionMatrices.gaussquad, Tuple{typeof(shape), Int}) ||
        throw(ArgumentError("shape must provide TransitionMatrices.gaussquad(shape, ngauss)"))
    hasproperty(shape, :m) ||
        throw(ArgumentError("shape must expose a relative refractive index field named `m`"))
    return nothing
end

function riccati_values(rho, nmax::Integer; hankel::Bool)
    psi = zeros(typeof(rho), nmax)
    psip = similar(psi)
    extra = TransitionMatrices.estimate_ricattibesselj_extra_terms(nmax, rho)
    z = zeros(typeof(rho), nmax + extra)
    TransitionMatrices.ricattibesselj!(psi, psip, z, nmax, extra, rho)

    if hankel
        chi = zeros(typeof(rho), nmax)
        chip = similar(chi)
        TransitionMatrices.ricattibessely!(chi, chip, nmax, rho)
        return CT.(psi .+ 1im .* chi), CT.(psip .+ 1im .* chip)
    end

    return CT.(psi), CT.(psip)
end

function angular_values(m::Integer, nmax::Integer, theta::Real)
    d, tau = TransitionMatrices.wigner_d_recursion(Float64, 0, m, nmax, theta;
                                                   deriv = true)
    pi_vals = zeros(Float64, nmax)
    for n in max(1, m):nmax
        pi_vals[n] = TransitionMatrices.pi_func(Float64, m, n, theta; d = d[n])
    end
    return d, pi_vals, tau
end

function mn_components(n::Integer, m::Integer, rho, psi_n, psip_n,
                       d_n, pi_n, tau_n)
    z = psi_n / rho
    dz = psip_n / rho

    m_r = zero(CT)
    m_theta = 1im * pi_n * z
    m_phi = -tau_n * z

    n_r = n * (n + 1) * d_n * z / rho
    n_theta = tau_n * dz
    n_phi = 1im * pi_n * dz

    return (M = (r = m_r, theta = m_theta, phi = m_phi),
            N = (r = n_r, theta = n_theta, phi = n_phi))
end

function add_mode_columns!(mat::Matrix{CT}, row0::Integer, col::Integer,
                           comps, qmedium, r::Real, rp::Real, sign::CT)
    M = comps.M
    N = comps.N

    # p = 1: electric M, magnetic q*N.
    mat[row0 + 1, col] += sign * (rp * M.r + r * M.theta)
    mat[row0 + 2, col] += sign * M.phi
    mat[row0 + 3, col] += sign * qmedium * (rp * N.r + r * N.theta)
    mat[row0 + 4, col] += sign * qmedium * N.phi

    # p = 2: electric N, magnetic q*M.
    mat[row0 + 1, col + 1] += sign * (rp * N.r + r * N.theta)
    mat[row0 + 2, col + 1] += sign * N.phi
    mat[row0 + 3, col + 1] += sign * qmedium * (rp * M.r + r * M.theta)
    mat[row0 + 4, col + 1] += sign * qmedium * M.phi

    return mat
end

function boundary_system_m(shape, lambda::Real, nmax::Integer,
                           m::Integer, ngauss::Integer)
    _check_axisymmetric_shape(shape)
    0 <= m <= nmax || throw(ArgumentError("m must satisfy 0 <= m <= nmax"))

    nmin = max(1, m)
    ns = collect(nmin:nmax)
    nn = length(ns)
    rows = 4 * ngauss
    unknowns = 4 * nn
    incident_cols = 2 * nn

    A = zeros(CT, rows, unknowns)
    B = zeros(CT, rows, incident_cols)
    weights = zeros(Float64, rows)

    xnodes, wnodes, radii, radius_derivatives = TransitionMatrices.gaussquad(shape, ngauss)
    k1 = 2pi / lambda
    q1 = CT(1)
    q2 = CT(getproperty(shape, :m))

    for ig in 1:ngauss
        theta = acos(Float64(xnodes[ig]))
        r = Float64(radii[ig])
        rp = Float64(radius_derivatives[ig])
        row0 = 4 * (ig - 1)

        surface_weight = Float64(wnodes[ig]) * r * sqrt(r^2 + rp^2)
        weights[(row0 + 1):(row0 + 4)] .= surface_weight

        d, pi_vals, tau = angular_values(m, nmax, theta)

        rho1 = k1 * r
        rho2 = getproperty(shape, :m) * k1 * r
        psi_j1, psip_j1 = riccati_values(rho1, nmax; hankel = false)
        psi_h1, psip_h1 = riccati_values(rho1, nmax; hankel = true)
        psi_j2, psip_j2 = riccati_values(rho2, nmax; hankel = false)

        for (local_idx, n) in enumerate(ns)
            inc_col = 2 * local_idx - 1
            sca_col = inc_col
            int_col = 2 * nn + inc_col

            inc_comps = mn_components(n, m, rho1, psi_j1[n], psip_j1[n],
                                      d[n], pi_vals[n], tau[n])
            sca_comps = mn_components(n, m, rho1, psi_h1[n], psip_h1[n],
                                      d[n], pi_vals[n], tau[n])
            int_comps = mn_components(n, m, rho2, psi_j2[n], psip_j2[n],
                                      d[n], pi_vals[n], tau[n])

            add_mode_columns!(B, row0, inc_col, inc_comps, q1, r, rp, CT(1))
            add_mode_columns!(A, row0, sca_col, sca_comps, q1, r, rp, CT(1))
            add_mode_columns!(A, row0, int_col, int_comps, q2, r, rp, CT(-1))
        end
    end

    return BoundarySystemM(A, B, weights, Int(m), nmin, Int(nmax))
end

function row_weight(system::BoundarySystemM)
    return sqrt.(system.weights ./ maximum(system.weights))
end

function weighted_blocks(system::BoundarySystemM)
    scale = row_weight(system)
    return system.A .* scale, system.B .* scale
end

function blocked_order(interleaved::AbstractMatrix{<:Number})
    nn2 = size(interleaved, 1)
    nn2 == size(interleaved, 2) ||
        throw(DimensionMismatch("mode block must be square"))
    iseven(nn2) || throw(DimensionMismatch("mode block size must be even"))
    nn = div(nn2, 2)
    out = zeros(CT, nn2, nn2)
    for ni in 1:nn, pi in 1:2
        ib = (pi - 1) * nn + ni
        ii = 2 * (ni - 1) + pi
        for nj in 1:nn, pj in 1:2
            jb = (pj - 1) * nn + nj
            ji = 2 * (nj - 1) + pj
            out[ib, jb] = interleaved[ii, ji]
        end
    end
    return out
end

function interleaved_order(blocked::AbstractMatrix{<:Number})
    nn2 = size(blocked, 1)
    nn2 == size(blocked, 2) ||
        throw(DimensionMismatch("mode block must be square"))
    iseven(nn2) || throw(DimensionMismatch("mode block size must be even"))
    nn = div(nn2, 2)
    out = zeros(CT, nn2, nn2)
    for ni in 1:nn, pi in 1:2
        ib = (pi - 1) * nn + ni
        ii = 2 * (ni - 1) + pi
        for nj in 1:nn, pj in 1:2
            jb = (pj - 1) * nn + nj
            ji = 2 * (nj - 1) + pj
            out[ii, ji] = blocked[ib, jb]
        end
    end
    return out
end

function block_from_tmatrix(Tmat, m::Integer, nmax::Integer)
    nmin = max(1, m)
    nn = nmax - nmin + 1
    block = zeros(CT, 2 * nn, 2 * nn)
    for pp in 1:2, p in 1:2, np in nmin:nmax, n in nmin:nmax
        i = (p - 1) * nn + n - nmin + 1
        j = (pp - 1) * nn + np - nmin + 1
        block[i, j] = CT(Tmat[m, n, m, np, p, pp])
    end
    return block
end

function axisymmetric_tmatrix_from_blocks(blocks::Vector{Matrix{CT}}, nmax::Integer)
    return AxisymmetricTransitionMatrix{CT, nmax, typeof(blocks), Float64}(blocks)
end

function truncate_tmatrix(Tmat, nmax::Integer)
    blocks = Matrix{CT}[]
    for m in 0:nmax
        push!(blocks, block_from_tmatrix(Tmat, m, nmax))
    end
    return axisymmetric_tmatrix_from_blocks(blocks, nmax)
end
