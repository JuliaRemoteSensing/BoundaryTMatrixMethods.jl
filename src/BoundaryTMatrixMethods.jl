module BoundaryTMatrixMethods

using LinearAlgebra
using TransitionMatrices

const CT = ComplexF64

include("boundary_system.jl")
include("solvers.jl")
include("diagnostics.jl")
include("correction.jl")
include("selector.jl")
include("metrics.jl")

export BoundarySystemM,
       ResidualDiagnostics,
       PmmApplicability,
       CorrectionResult,
       SelectionResult,
       boundary_system_m,
       svm_tmatrix,
       pmm_tmatrix,
       solve_svm_block,
       solve_pmm_block,
       boundary_residual,
       projected_boundary_residual,
       residual_diagnostics,
       diagnose_pmm_applicability,
       residual_corrected_iitm,
       pmm_residual_corrected_tmatrix,
       pmm_guided_iitm,
       auto_tmatrix,
       trust_blend_tmatrix,
       truncate_tmatrix,
       relative_tmatrix_error,
       mueller_relative_error,
       cross_section_relative_errors

end
