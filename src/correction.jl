struct CorrectionResult
    T::Any
    base::Any
    proposal::Any
    accepted::Bool
    method::Symbol
    diagnostics::PmmApplicability
    train_residual_before::Float64
    train_residual_after::Float64
    relative_step::Float64
    reasons::Vector{String}
end

function strongest_rows(column::AbstractVector{<:Number}, max_rows::Integer)
    max_rows >= length(column) && return collect(eachindex(column))
    order = sortperm(abs2.(column); rev = true)
    selected = order[1:max_rows]
    sort!(selected)
    return selected
end

function pmm_residual_correct_block(system::BoundarySystemM,
                                    base_blocked::AbstractMatrix{<:Number};
                                    max_rows_per_column::Integer = 8,
                                    regularization::Real = 1e-4)
    _, B, C, D, nn2 = _split_boundary_blocks(system)
    base = interleaved_order(base_blocked)
    corrected = copy(base)
    before2 = 0.0
    after2 = 0.0
    project = _external_projector(D)

    for col in 1:nn2
        t0 = view(base, :, col)
        b = view(B, :, col)
        r0 = project(C * t0 .+ b)
        before2 += abs2(norm(r0))

        selected = strongest_rows(t0, max_rows_per_column)
        G = project(C[:, selected])
        lhs = G' * G + CT(regularization) * Matrix{CT}(I, length(selected),
                                                       length(selected))
        rhs = -(G' * r0)
        delta = lhs \ rhs
        corrected[selected, col] .+= delta

        t1 = view(corrected, :, col)
        r1 = project(C * t1 .+ b)
        after2 += abs2(norm(r1))
    end

    return (block = blocked_order(corrected),
            train_residual_before = sqrt(before2) / max(norm(B), eps(Float64)),
            train_residual_after = sqrt(after2) / max(norm(B), eps(Float64)))
end

function pmm_residual_corrected_tmatrix(base, shape, lambda::Real, nmax::Integer,
                                        ngauss::Integer;
                                        max_rows_per_column::Integer = 8,
                                        regularization::Real = 1e-4)
    blocks = Matrix{CT}[]
    before2 = 0.0
    after2 = 0.0
    for m in 0:nmax
        system = boundary_system_m(shape, lambda, nmax, m, ngauss)
        base_block = block_from_tmatrix(base, m, nmax)
        corrected = pmm_residual_correct_block(system, base_block;
                                               max_rows_per_column = max_rows_per_column,
                                               regularization = regularization)
        push!(blocks, corrected.block)
        before2 += corrected.train_residual_before^2
        after2 += corrected.train_residual_after^2
    end

    return (T = axisymmetric_tmatrix_from_blocks(blocks, nmax),
            train_residual_before = sqrt(before2),
            train_residual_after = sqrt(after2))
end

function trust_blend_tmatrix(left, right, alpha::Real; nmax::Integer = min(size(left, 2),
                                                                          size(right, 2)))
    0 <= alpha <= 1 || throw(ArgumentError("alpha must lie in [0, 1]"))
    blocks = Matrix{CT}[]
    for m in 0:nmax
        L = block_from_tmatrix(left, m, nmax)
        R = block_from_tmatrix(right, m, nmax)
        push!(blocks, (1 - alpha) .* L .+ alpha .* R)
    end
    return axisymmetric_tmatrix_from_blocks(blocks, nmax)
end

function residual_corrected_iitm(base, shape, lambda::Real, nmax::Integer;
                                 ngauss::Integer = 120,
                                 max_rows_per_column::Integer = 8,
                                 regularization::Real = 1e-4,
                                 gate::Union{Nothing, PmmApplicability} = nothing,
                                 step_tol::Real = 0.75,
                                 gate_kwargs...)
    diagnostics = isnothing(gate) ?
        diagnose_pmm_applicability(shape, lambda, nmax; ngauss = ngauss, gate_kwargs...) :
        gate

    proposal = pmm_residual_corrected_tmatrix(base, shape, lambda, nmax, ngauss;
                                              max_rows_per_column = max_rows_per_column,
                                              regularization = regularization)
    step = relative_tmatrix_error(proposal.T, base; nmax = nmax)
    reasons = copy(diagnostics.reasons)
    step <= step_tol ||
        push!(reasons, "PMM-guided correction step is larger than the trust radius")

    accepted = diagnostics.accepted && step <= step_tol
    Tfinal = accepted ? proposal.T : base
    method = accepted ? :pmm_residual_corrected_iitm : :iitm_base
    return CorrectionResult(Tfinal, base, proposal.T, accepted, method, diagnostics,
                            Float64(proposal.train_residual_before),
                            Float64(proposal.train_residual_after),
                            Float64(step), reasons)
end

function pmm_guided_iitm(shape, lambda::Real, nmax::Integer;
                         Nr::Integer = 6,
                         Ntheta::Integer = 40,
                         ngauss::Integer = 120,
                         kwargs...)
    base = transition_matrix_iitm(shape, lambda, nmax, Nr, Ntheta)
    return residual_corrected_iitm(base, shape, lambda, nmax;
                                   ngauss = ngauss, kwargs...)
end
