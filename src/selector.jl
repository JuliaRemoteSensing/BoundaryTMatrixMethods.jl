struct SelectionResult
    T::Any
    method::Symbol
    diagnostics::PmmApplicability
    fallback::Any
    reasons::Vector{String}
end

function auto_tmatrix(shape, lambda::Real, nmax::Integer;
                      Nr::Integer = 6,
                      Ntheta::Integer = 40,
                      pmm_ngauss::Integer = 120,
                      prefer::Symbol = :pmm,
                      kwargs...)
    prefer in (:pmm, :iitm) ||
        throw(ArgumentError("prefer must be either :pmm or :iitm"))

    fallback = transition_matrix_iitm(shape, lambda, nmax, Nr, Ntheta)
    diagnostics = diagnose_pmm_applicability(shape, lambda, nmax;
                                             ngauss = pmm_ngauss, kwargs...)

    if prefer == :pmm && diagnostics.accepted
        return SelectionResult(diagnostics.T, :pmm, diagnostics, fallback, String[])
    end

    reasons = diagnostics.accepted ? ["user preference selected IITM fallback"] :
              copy(diagnostics.reasons)
    return SelectionResult(fallback, :iitm, diagnostics, fallback, reasons)
end
