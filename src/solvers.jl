function solve_svm_block(system::BoundarySystemM; row_weighted::Bool = true,
                         regularization::Real = 0.0)
    A = system.A
    B = system.B
    if row_weighted
        scale = row_weight(system)
        A = A .* scale
        B = B .* scale
    end

    if regularization == 0
        X = -(A \ B)
    else
        lhs = A' * A + CT(regularization) * Matrix{CT}(I, size(A, 2), size(A, 2))
        rhs = A' * B
        X = -(lhs \ rhs)
    end

    nn2 = size(B, 2)
    return blocked_order(X[1:nn2, :])
end

function solve_pmm_block(system::BoundarySystemM; regularization::Real = 0.0)
    A, B = weighted_blocks(system)
    lhs = A' * A
    rhs = A' * B
    if regularization != 0
        lhs += CT(regularization) * Matrix{CT}(I, size(lhs, 1), size(lhs, 2))
    end
    X = -(lhs \ rhs)

    nn2 = size(B, 2)
    return blocked_order(X[1:nn2, :])
end

function _tmatrix_from_boundary_blocks(shape, lambda::Real, nmax::Integer,
                                       ngauss::Integer, block_solver; kwargs...)
    blocks = Matrix{CT}[]
    for m in 0:nmax
        system = boundary_system_m(shape, lambda, nmax, m, ngauss)
        push!(blocks, block_solver(system; kwargs...))
    end
    return axisymmetric_tmatrix_from_blocks(blocks, nmax)
end

function svm_tmatrix(shape, lambda::Real, nmax::Integer, ngauss::Integer;
                     row_weighted::Bool = true,
                     regularization::Real = 0.0)
    return _tmatrix_from_boundary_blocks(shape, lambda, nmax, ngauss,
                                         solve_svm_block;
                                         row_weighted = row_weighted,
                                         regularization = regularization)
end

function pmm_tmatrix(shape, lambda::Real, nmax::Integer, ngauss::Integer;
                     regularization::Real = 0.0)
    return _tmatrix_from_boundary_blocks(shape, lambda, nmax, ngauss,
                                         solve_pmm_block;
                                         regularization = regularization)
end
