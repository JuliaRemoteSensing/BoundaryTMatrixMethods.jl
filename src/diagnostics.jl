struct ResidualDiagnostics
    residual::Float64
    projected_residual::Float64
    csca::Float64
    cext::Float64
    absorption::Float64
end

struct PmmApplicability
    accepted::Bool
    status::Symbol
    reasons::Vector{String}
    nmax::Int
    ngauss::Int
    refined_nmax::Int
    refined_ngauss::Int
    convergence_error::Float64
    residual::Float64
    projected_residual::Float64
    energy_violation::Float64
    csca::Float64
    cext::Float64
    T::Any
end

function _split_boundary_blocks(system::BoundarySystemM)
    A, B = weighted_blocks(system)
    nn2 = size(B, 2)
    C = view(A, :, 1:nn2)
    D = view(A, :, (nn2 + 1):(2 * nn2))
    return A, B, C, D, nn2
end

function _external_projector(D)
    qmat = Matrix(qr(Matrix(D)).Q)[:, 1:size(D, 2)]
    return v -> v .- qmat * (qmat' * v)
end

function boundary_residual(system::BoundarySystemM, blocked_sca::AbstractMatrix{<:Number})
    _, B, C, D, nn2 = _split_boundary_blocks(system)
    X = interleaved_order(blocked_sca)
    rhs = -(C * X .+ B)
    Y = D \ rhs
    residual = C * X .+ D * Y .+ B
    return norm(residual) / max(norm(B), eps(Float64))
end

function projected_boundary_residual(system::BoundarySystemM,
                                     blocked_sca::AbstractMatrix{<:Number})
    _, B, C, D, _ = _split_boundary_blocks(system)
    X = interleaved_order(blocked_sca)
    project = _external_projector(D)
    total2 = 0.0
    for col in axes(B, 2)
        r = project(C * view(X, :, col) .+ view(B, :, col))
        total2 += abs2(norm(r))
    end
    return sqrt(total2) / max(norm(B), eps(Float64))
end

function boundary_residual(Tmat, shape, lambda::Real, nmax::Integer,
                           ngauss::Integer)
    total2 = 0.0
    for m in 0:nmax
        system = boundary_system_m(shape, lambda, nmax, m, ngauss)
        block = block_from_tmatrix(Tmat, m, nmax)
        total2 += boundary_residual(system, block)^2
    end
    return sqrt(total2)
end

function projected_boundary_residual(Tmat, shape, lambda::Real, nmax::Integer,
                                     ngauss::Integer)
    total2 = 0.0
    for m in 0:nmax
        system = boundary_system_m(shape, lambda, nmax, m, ngauss)
        block = block_from_tmatrix(Tmat, m, nmax)
        total2 += projected_boundary_residual(system, block)^2
    end
    return sqrt(total2)
end

function residual_diagnostics(Tmat, shape, lambda::Real, nmax::Integer,
                              ngauss::Integer)
    csca = Float64(scattering_cross_section(Tmat, lambda))
    cext = Float64(extinction_cross_section(Tmat, lambda))
    return ResidualDiagnostics(
        boundary_residual(Tmat, shape, lambda, nmax, ngauss),
        projected_boundary_residual(Tmat, shape, lambda, nmax, ngauss),
        csca,
        cext,
        cext - csca,
    )
end

function relative_tmatrix_error(candidate, reference; nmax::Integer = min(size(candidate, 2),
                                                                         size(reference, 2)))
    num2 = 0.0
    den2 = 0.0
    for m in 0:nmax
        A = block_from_tmatrix(candidate, m, nmax)
        B = block_from_tmatrix(reference, m, nmax)
        num2 += sum(abs2, A .- B)
        den2 += sum(abs2, B)
    end
    return sqrt(num2) / max(sqrt(den2), eps(Float64))
end

function _energy_violation(Tmat, lambda::Real)
    csca = Float64(scattering_cross_section(Tmat, lambda))
    cext = Float64(extinction_cross_section(Tmat, lambda))
    return max(0.0, csca - cext) / max(abs(cext), eps(Float64))
end

function diagnose_pmm_applicability(shape, lambda::Real, nmax::Integer;
                                    ngauss::Integer = 120,
                                    refined_nmax::Integer = nmax + 2,
                                    refined_ngauss::Integer = max(ngauss + 40, ceil(Int, 1.5 * ngauss)),
                                    regularization::Real = 0.0,
                                    convergence_tol::Real = 5e-3,
                                    residual_tol::Real = 0.5,
                                    projected_residual_tol::Real = 0.5,
                                    energy_tol::Real = 1e-3)
    Tpmm = pmm_tmatrix(shape, lambda, nmax, ngauss;
                       regularization = regularization)
    Trefined = pmm_tmatrix(shape, lambda, refined_nmax, refined_ngauss;
                           regularization = regularization)
    Trefined_trunc = truncate_tmatrix(Trefined, nmax)
    convergence = relative_tmatrix_error(Tpmm, Trefined_trunc; nmax = nmax)
    diag = residual_diagnostics(Tpmm, shape, lambda, nmax, ngauss)
    energy = _energy_violation(Tpmm, lambda)

    reasons = String[]
    convergence <= convergence_tol ||
        push!(reasons, "PMM T-matrix did not converge under refined nmax/Ng")
    diag.residual <= residual_tol ||
        push!(reasons, "boundary residual is above the requested tolerance")
    diag.projected_residual <= projected_residual_tol ||
        push!(reasons, "projected boundary residual is above the requested tolerance")
    energy <= energy_tol ||
        push!(reasons, "cross sections violate the passive-particle energy gate")

    accepted = isempty(reasons)
    status = accepted ? :accepted : :rejected
    return PmmApplicability(accepted, status, reasons, Int(nmax), Int(ngauss),
                            Int(refined_nmax), Int(refined_ngauss),
                            Float64(convergence), diag.residual,
                            diag.projected_residual, Float64(energy),
                            diag.csca, diag.cext, Tpmm)
end
